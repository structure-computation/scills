#include <string>
#include <algorithm>

#include "mpi_lmt_functions.h"
#include "fonctions_affichage_interactivite.h"

using namespace LMT;
using namespace std;

///Update le param_comp->jeu pour les contact_jeu et les jeu_impose
template<class TV2>
void update_jeu(TV2 &Inter,Param &process){
  typedef typename TV2::T T;
  std::vector<Ex> symbols;
  if (TV2::dim==2) {
    symbols.push_back("x");
    symbols.push_back("y");
  } else if (TV2::dim==3) {
    symbols.push_back("x");
    symbols.push_back("y");
    symbols.push_back("z");
  }
  //if (process.rank==0) std::cout << Inter.param_comp->jeu << std::endl;
  Vec<string> jeu_cut=tokenize(Inter.param_comp->fcts_spatiales,';');
  if (jeu_cut.size() == 1) {//Jeu normal
    Ex expr;
    expr = read_ex(Inter.param_comp->fcts_spatiales.c_str(),symbols);
    for(unsigned k=0;k<Inter.side[0].nodeeq.size();++k) {
      T data;
      Ex::MapExNum var;
      for(unsigned d2=0;d2<TV2::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
        var[symbols[d2]]= Inter.side[0].nodeeq[k][d2];
      data = (T)expr.subs_numerical(var);
      Inter.param_comp->jeu[range(TV2::dim*k,TV2::dim*(k+1))]=data*Inter.side[0].neq[range(TV2::dim*k,TV2::dim*(k+1))];
    }
  } else {//Jeu complet
    Vec<Ex> expr;
    expr.resize(TV2::dim);

    for(unsigned d2=0;d2<TV2::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
      expr[d2] = read_ex(jeu_cut[d2],symbols);
    }
    for(unsigned k=0;k<Inter.side[0].nodeeq.size();++k) {
      Vec<T,TV2::dim> data;
      Ex::MapExNum var;
      for(unsigned d2=0;d2<TV2::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
        var[symbols[d2]]= Inter.side[0].nodeeq[k][d2];
      }
      for(unsigned d2=0;d2<TV2::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
        data[d2] = (T)expr[d2].subs_numerical(var);
      Inter.param_comp->jeu[range(TV2::dim*k,TV2::dim*(k+1))]=data;

    }
  }
  //if (process.rank==0) std::cout << Inter.param_comp->jeu << std::endl;
}

///Convertit une chaine de caractere en minuscule
inline void minuscule(string &s) {
    transform(s.begin(), s.end(), s.begin(), static_cast<int (*)(int)>(tolower));
}

inline Vec<string> get_token(istream &is) {
    string entree;
    getline(is,entree);
    Vec<string> vec_entree;
    vec_entree=tokenize(entree,' ');
    return vec_entree;
}

/// Nom interactif de fonction
typedef void (*parametre_fct)(Vec<string> &v,Param &process);

/// Initialisation des mots cles dans la map
inline void initialize_mots_cles(map<string, parametre_fct> &param_map) {
    param_map[ "evol" ] = evol;
    param_map[ "evol_inter" ] = evol_inter;
    //param_map[ "trac_sst" ] = trac_sst;
    //param_map[ "trac_inter" ] = trac_inter;
    param_map[ "trac_sst_temps" ] = trac_sst_temps;
    param_map[ "trac_inter_temps" ] = trac_inter_temps;
    param_map[ "mesh_inter" ] = trac_mesh_inter;
    param_map[ "mesh_sst" ] = trac_mesh_sst;
    param_map[ "trac_error" ] = trac_error;
    param_map[ "trac_ener" ] = trac_ener;
    param_map[ "trac_cl" ] = trac_cl;
    param_map[ "calcul" ] = calcul;
    param_map[ "modif_param" ] = modif_param;
    param_map[ "modif_param_inter" ] = modif_param_inter;
    param_map[ "help" ] = help;
}

/** \ingroup Post_Traitement
 \brief Fonction principale permettant une interactivite sur les résultats. 
 
 Si le mot clé interactivite est mis à 1 dans le fichier xml, on fait appel à cette fonction. On a alors accès à différents mots clés pour 
 - tracer la déformée en sélectionnant les champs souhaités
 - tracer le déplacement d'un point donné
 - tracer les vitesses, efforts ... sur les interfaces
 - tracer la courbe d'erreur latin
 - effectuer des opérations définis par l'utilisateur (exemple calcul d'énergie par interface)
 
 Pour accéder aux différents mots clés il suffit de taper help puis entrer.
 */
template<class TV1, class TV2,class TV3, class TV4,class GLOB, class TV5, class TV6>
void interactivite(TV1 &S,TV3 &SubS, TV2 &Inter, TV4 &SubI, Param &process, GLOB &Global, TV5 &CL, TV6 &n) {
    string str;
    LMT::Vec<string> v;

    v.resize(1);
    v[0]="debut";
    static map<string, parametre_fct> param_map;

    // initialiser la map si ce n'est pas fait
    if ( param_map.empty() )
        initialize_mots_cles(param_map);

    std::ifstream f(process.affichage->command_file.data());

    string entree;
    while( (v.size()==0 or v[0] != "exit") ) {
        if (process.size >1){
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(1);
        }
        if (process.rank == 0)
            std::cout << "> ";
        //std::cout << process.rank << " : "  << process.size << std::endl;

        if (process.rank == 0) {
            if (f) {
                getline(f,str);
                std::cout << str << std::endl;
            } else {
                getline(std::cin,str);
                process.affichage->command_file="";///permet de permettre le reaffichage en mode interactivite a la fin de la lecture du fichier de commande
            }
        }
        if (process.size > 1) MPI_Bcast(str,0);

        // rechercher la fonction associée à Param
        v=tokenize(str,' ');
        //std::cout << process.rank << " : " << v << std::endl;
        if( v.size() > 0 and v[0] != "exit") {
            map<string, parametre_fct>::const_iterator i = param_map.find( v[0] );
            if ( i == param_map.end() ) {
                // échec
                if (process.rank == 0)
                    parametre_inconnu(v,process);
            } else {
                minuscule(v[0]);
                // appeler la fonction associée pour regler les parametres
                (i->second)(v,process);
                //appeler la fonction a executer
                if(v[0]=="evol")
                    affichage_depl_pt(SubS, process);
                if(v[0]=="evol_inter")
                    affichage_var_inter(SubS,Inter, process);
                //if(v[0]=="trac_sst") {affichage_resultats(SubS,process); process.affichage->fichiers_paraview_sst_crees=1;}
                if(v[0]=="trac_sst_temps") {
                    if (process.latin->save_depl_SST!=1){
                      if(process.rank == 0) std::cout << "ATTENTION il faut mettre save_depl_SST a 1 pour utiliser cette commande" << std::endl;
                    } else {
                    if(process.affichage->fichiers_paraview_sst_crees==1) {
                        if (process.rank == 0)
                            affichage_resultats_temps(process);
                    } else {
                        affichage_resultats(SubS,process);
                        if (process.size >1)
                            MPI_Barrier(MPI_COMM_WORLD);
                        process.affichage->fichiers_paraview_sst_crees=1;
                        if (process.rank == 0)
                            affichage_resultats_temps(process);
                    }
                    }
                }
                //if (process.rank == 0) if(v[0]=="trac_inter") affichage_inter_data(Inter, S, process);
                if(v[0]=="trac_inter_temps") {
                  if(process.affichage->fichiers_paraview_inter_crees==1){
                        if (process.rank == 0)
                          affichage_inter_temps(process);
                  } else {
                            affichage_resultats_inter(SubI,S,process);
                            if (process.size >1)
                                MPI_Barrier(MPI_COMM_WORLD);
                            process.affichage->fichiers_paraview_inter_crees=1;
                            if (process.rank == 0)
                                affichage_inter_temps(process);
                        }
                }
                if(v[0]=="mesh_sst") affichage_maillage(SubS,SubI,S,process);
                if(v[0]=="mesh_inter") affichage_maillage(SubS,SubI,S,process);
                if (process.rank == 0)
                    if(v[0]=="trac_error") {
                        std::cout << process.latin->error[range(process.latin->iter+1)] << std::endl;
                        if (process.affichage->command_file=="") {
                            GnuPlot gp;
                            gp.plot(log10(process.latin->error[range(process.latin->iter+1)]));
                            gp.wait();
                        }
                    }
                if(v[0]=="trac_ener")
                    affichage_energie(SubS,Inter,process);
                if(v[0]=="trac_cl")
                  affichage_cl(SubS,Inter,process); if (process.size >1) MPI_Barrier(MPI_COMM_WORLD);
                if(v[0]=="modif_param_inter"){
                  if (v.size()> 1) {
                  if (find(range(Inter.size()),LMT::_1==(unsigned)atoi(v[1].c_str()))){
                    unsigned q=atoi(v[1].c_str());
                    if (Inter[q].comp=="effort" or Inter[q].comp=="depl" or Inter[q].comp=="depl_normal"){
                      if (process.rank==0) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction spatiale remplacant " << CL[Inter[q].refCL].fcts_spatiales << std::endl;
                      if (process.rank == 0) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.size > 1) MPI_Bcast(str,0);
                      Vec<string> fctlu=tokenize(str,';');
                      if (fctlu.size()==CL[Inter[q].refCL].fcts_spatiales.size()) CL[Inter[q].refCL].fcts_spatiales=fctlu;
                      else if (process.rank==0) std::cout << "La fonction donnée n'est pas sous la bonne forme : 0;0;0" << std::endl;
                      
                      if (v.size()==3 and v[2]=="temps"){
                          for( unsigned kk=0;kk<CL[Inter[q].refCL].fcts_temporelles.size() ;kk++ ){
                              if (process.rank==0) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction temporelle remplacant " << CL[Inter[q].refCL].fcts_temporelles[kk] << std::endl;
                              if (process.rank == 0) {
                                  if (f) {
                                      getline(f,str);
                                      std::cout << str << std::endl;
                                  } else {
                                      getline(std::cin,str);
                                      process.affichage->command_file="";
                                  }
                              }
                              if (process.size > 1) MPI_Bcast(str,0);
                              CL[Inter[q].refCL].fcts_temporelles[kk]=tokenize(str,';');
                          //Vec<string> fctlu=tokenize(str,';');
                          //if (fctlu.size()==CL[Inter[q].refCL].fcts_spatiales.size()) CL[Inter[q].refCL].fcts_spatiales=fctlu;
                          //else if (process.rank==0) std::cout << "La fonction donnée n'est pas sous la bonne forme : 0;0;0" << std::endl;
                          }
                      }
                    } else if (Inter[q].comp=="Jeu_impose" or Inter[q].comp=="Contact_jeu"){
                      if (process.rank==0) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction spatiale remplacant " << Inter[q].param_comp->fcts_spatiales << std::endl;
                      if (process.rank == 0) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.size > 1) MPI_Bcast(str,0);
                      Vec<string> t1,t2;
                      t1=tokenize(str,';');
                      t2=tokenize(Inter[q].param_comp->fcts_spatiales,';');
                      if (t1.size()== t2.size() or t1.size()==1){
                      Inter[q].param_comp->fcts_spatiales=str;
                      if (find(process.multi_mpi->listinter,LMT::_1==q)) update_jeu(Inter[q],process);
                      } else if (process.rank==0) std::cout << "La fonction donnée n'est pas sous la bonne forme : 0;0;0" << std::endl;
                    } else if (Inter[q].comp=="Contact" or Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique"){
                      if (process.rank==0) std::cout << "Interface de type " << Inter[q].comp << " donner le nouveau coefficient de frottement " << Inter[q].param_comp->coeffrottement << std::endl;
                      if (process.rank == 0) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.size > 1) MPI_Bcast(str,0);
                      Inter[q].param_comp->coeffrottement=atof(str.c_str());
                    } else {
                      if (process.rank==0) std::cout << "Interface de type " << Inter[q].comp << " non modifiable" << std::endl;
                    }
                  } else {
                    if (process.rank==0) std::cout << "mauvais numero d interface" << std::endl;
                  }
                  } else {
                    if (process.rank==0) std::cout << "donner un numero d interface" << std::endl;
                  }
                }
                if(v[0]=="calcul") {
                    process.affichage->fichiers_paraview_inter_crees=0;
                    process.affichage->fichiers_paraview_sst_crees=0;
                    process.reprise_calcul=2;
                  if(process.nom_calcul=="incr") multiscale_iterate_incr(S,SubS, Inter,SubI, process, Global, CL,n);
                  else if(process.nom_calcul=="latin") multiscale_iterate_latin(S,SubS, Inter,SubI, process, Global, CL);

                }
            }
        }
    }
}


