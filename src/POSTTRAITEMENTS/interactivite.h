#include "../UTILS/Sc2String.h"
#include <algorithm>

#include "ParallelisationData.h"
#include "mpi_lmt_functions.h"
#include "fonctions_affichage_interactivite.h"

using namespace LMT;


///Update le param_comp->jeu pour les contact_jeu et les jeu_impose
template<class TV2>
void update_jeu(TV2 &Inter,Process &process){
    std::vector<Ex> symbols;
    if (DIM==2) {
        symbols.push_back("x");
        symbols.push_back("y");
    } else if (DIM==3) {
        symbols.push_back("x");
        symbols.push_back("y");
        symbols.push_back("z");
    }
    //if (process.parallelisation->is_master_cpu()) std::cout << Inter.param_comp->jeu << std::endl;
    Vec<Sc2String> jeu_cut=tokenize(Inter.param_comp->fcts_spatiales,';');
    if (jeu_cut.size() == 1) {//Jeu normal
        Ex expr;
        expr = read_ex(Inter.param_comp->fcts_spatiales.c_str(),symbols);
        for(unsigned k=0;k<Inter.side[0].nodeeq.size();++k) {
            TYPEREEL data;
            Ex::MapExNum var;
            for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= Inter.side[0].nodeeq[k][d2];
            data = (TYPEREEL)expr.subs_numerical(var);
            Inter.param_comp->jeu[range(DIM*k,DIM*(k+1))]=data*Inter.side[0].neq[range(DIM*k,DIM*(k+1))];
        }
    } else {//Jeu complet
        Vec<Ex> expr;
        expr.resize(DIM);

        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(jeu_cut[d2],symbols);
        }
        for(unsigned k=0;k<Inter.side[0].nodeeq.size();++k) {
            Vec<TYPEREEL,DIM> data;
            Ex::MapExNum var;
            for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= Inter.side[0].nodeeq[k][d2];
            for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                data[d2] = (TYPEREEL)expr[d2].subs_numerical(var);
            Inter.param_comp->jeu[range(DIM*k,DIM*(k+1))]=data;

        }
    }
    //if (process.parallelisation->is_master_cpu()) std::cout << Inter.param_comp->jeu << std::endl;
}

///Convertit une chaine de caractere en minuscule
inline void minuscule(Sc2String &s) {
    transform(s.begin(), s.end(), s.begin(), static_cast<int (*)(int)>(tolower));
}

inline Vec<Sc2String> get_token(istream &is) {
    Sc2String entree;
    getline(is,entree);
    Vec<Sc2String> vec_entree;
    vec_entree=tokenize(entree,' ');
    return vec_entree;
}

/// Nom interactif de fonction
typedef void (*parametre_fct)(Vec<Sc2String> &v,Process &process);

/// Initialisation des mots cles dans la map
inline void initialize_mots_cles(map<Sc2String, parametre_fct> &param_map) {
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
 \brief Fonction principale permettant une interactivite sur les r�sultats. 
 
 Si le mot cl� interactivite est mis � 1 dans le fichier xml, on fait appel � cette fonction. On a alors acc�s � diff�rents mots cl�s pour 
 - tracer la d�form�e en s�lectionnant les champs souhait�s
 - tracer le d�placement d'un point donn�
 - tracer les vitesses, efforts ... sur les interfaces
 - tracer la courbe d'erreur latin
 - effectuer des op�rations d�finis par l'utilisateur (exemple calcul d'�nergie par interface)
 
 Pour acc�der aux diff�rents mots cl�s il suffit de taper help puis entrer.
 */
template<class TV1, class TV2,class TV3, class TV4,class GLOB, class TV5, class TV6>
void interactivite(TV1 &S,TV3 &SubS, TV2 &Inter, TV4 &SubI, Process &process, GLOB &Global, TV5 &CL, TV6 &n) {
    Sc2String str;
    LMT::Vec<Sc2String> v;

    v.resize(1);
    v[0]="debut";
    static map<Sc2String, parametre_fct> param_map;

    // initialiser la map si ce n'est pas fait
    if ( param_map.empty() )
        initialize_mots_cles(param_map);

    std::ifstream f(process.affichage->command_file.data());

    Sc2String entree;
    while( (v.size()==0 or v[0] != "exit") ) {
        process.parallelisation->synchronisation();
        if(process.parallelisation->is_multi_cpu())
            sleep(1);
        }
        if (process.parallelisation->is_master_cpu())
            std::cout << "> ";
        //std::cout << process.parallelisation->rank << " : "  << process.parallelisation->size << std::endl;

        if (process.parallelisation->is_master_cpu()) {
            if (f) {
                getline(f,str);
                std::cout << str << std::endl;
            } else {
                getline(std::cin,str);
                process.affichage->command_file="";///permet de permettre le reaffichage en mode interactivite a la fin de la lecture du fichier de commande
            }
        }
        if (process.parallelisation->is_multi_cpu()) MPI_Bcast(str,0);

        // rechercher la fonction associ�e � Process
        v=tokenize(str,' ');
        //std::cout << process.parallelisation->rank << " : " << v << std::endl;
        if( v.size() > 0 and v[0] != "exit") {
            map<Sc2String, parametre_fct>::const_iterator i = param_map.find( v[0] );
            if ( i == param_map.end() ) {
                // �chec
                if (process.parallelisation->is_master_cpu())
                    parametre_inconnu(v,process);
            } else {
                minuscule(v[0]);
                // appeler la fonction associ�e pour regler les parametres
                (i->second)(v,process);
                //appeler la fonction a executer
                if(v[0]=="evol")
                    affichage_depl_pt(SubS, process);
                if(v[0]=="evol_inter")
                    affichage_var_inter(SubS,Inter, process);
                //if(v[0]=="trac_sst") {affichage_resultats(SubS,process); process.affichage->fichiers_paraview_sst_crees=1;}
                if(v[0]=="trac_sst_temps") {
                    if (process.latin->save_depl_SST!=1){
                      if(process.parallelisation->is_master_cpu()) std::cout << "ATTENTION il faut mettre save_depl_SST a 1 pour utiliser cette commande" << std::endl;
                    } else {
                    if(process.affichage->fichiers_paraview_sst_crees==1) {
                        if (process.parallelisation->is_master_cpu())
                            affichage_resultats_temps(process);
                    } else {
                        affichage_resultats(SubS,process);
                        process.parallelisation->synchronisation();
                        process.affichage->fichiers_paraview_sst_crees=1;
                        if (process.parallelisation->is_master_cpu())
                            affichage_resultats_temps(process);
                    }
                    }
                }
                //if (process.parallelisation->is_master_cpu()) if(v[0]=="trac_inter") affichage_inter_data(Inter, S, process);
                if(v[0]=="trac_inter_temps") {
                  if(process.affichage->fichiers_paraview_inter_crees==1){
                        if (process.parallelisation->is_master_cpu())
                          affichage_inter_temps(process);
                  } else {
                            affichage_resultats_inter(SubI,S,process);
                            process.parallelisation->synchronisation();
                            process.affichage->fichiers_paraview_inter_crees=1;
                            if (process.parallelisation->is_master_cpu())
                                affichage_inter_temps(process);
                        }
                }
                if(v[0]=="mesh_sst") affichage_maillage(SubS,SubI,S,process);
                if(v[0]=="mesh_inter") affichage_maillage(SubS,SubI,S,process);
                if (process.parallelisation->is_master_cpu())
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
                  affichage_cl(SubS,Inter,process);
                process.parallelisation->synchronisation();
                if(v[0]=="modif_param_inter"){
                  if (v.size()> 1) {
                  if (find(range(Inter.size()),LMT::_1==(unsigned)atoi(v[1].c_str()))){
                    unsigned q=atoi(v[1].c_str());
                    if (Inter[q].comp=="depl" or Inter[q].comp=="depl_normal" or Inter[q].comp=="depl_nul" or Inter[q].comp=="vit_nulle" or Inter[q].comp=="vit" or Inter[q].comp=="vit_normale" or Inter[q].comp=="effort" or Inter[q].comp=="effort_normal"){
                      if (process.parallelisation->is_master_cpu()) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction spatiale remplacant " << CL[Inter[q].refCL].fcts_spatiales << std::endl;
                      if (process.parallelisation->is_master_cpu()) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.parallelisation->is_multi_cpu()) MPI_Bcast(str,0);
                      Vec<Sc2String> fctlu=tokenize(str,';');
                      if (fctlu.size()==CL[Inter[q].refCL].fcts_spatiales.size()) CL[Inter[q].refCL].fcts_spatiales=fctlu;
                      else if (process.parallelisation->is_master_cpu()) std::cout << "La fonction donn�e n'est pas sous la bonne forme : 0;0;0" << std::endl;
                      
                      if (v.size()==3 and v[2]=="temps"){
                          for( unsigned kk=0;kk<CL[Inter[q].refCL].fcts_temporelles.size() ;kk++ ){
                              if (process.parallelisation->is_master_cpu()) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction temporelle remplacant " << CL[Inter[q].refCL].fcts_temporelles[kk] << std::endl;
                              if (process.parallelisation->is_master_cpu()) {
                                  if (f) {
                                      getline(f,str);
                                      std::cout << str << std::endl;
                                  } else {
                                      getline(std::cin,str);
                                      process.affichage->command_file="";
                                  }
                              }
                              if (process.parallelisation->is_multi_cpu()) MPI_Bcast(str,0);
                              CL[Inter[q].refCL].fcts_temporelles[kk]=str;
                          //Vec<Sc2String> fctlu=tokenize(str,';');
                          //if (fctlu.size()==CL[Inter[q].refCL].fcts_spatiales.size()) CL[Inter[q].refCL].fcts_spatiales=fctlu;
                          //else if (process.parallelisation->is_master_cpu()) std::cout << "La fonction donn�e n'est pas sous la bonne forme : 0;0;0" << std::endl;
                          }
                      }
                    } else if (Inter[q].comp=="Jeu_impose" or Inter[q].comp=="Contact_jeu"){
                      if (process.parallelisation->is_master_cpu()) std::cout << "Interface de type " << Inter[q].comp << " donner la nouvelle fonction spatiale remplacant " << Inter[q].param_comp->fcts_spatiales << std::endl;
                      if (process.parallelisation->is_master_cpu()) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.parallelisation->is_multi_cpu()) MPI_Bcast(str,0);
                      Vec<Sc2String> t1,t2;
                      t1=tokenize(str,';');
                      t2=tokenize(Inter[q].param_comp->fcts_spatiales,';');
                      if (t1.size()== t2.size() or t1.size()==1){
                      Inter[q].param_comp->fcts_spatiales=str;
                      if (find(process.parallelisation->listinter,LMT::_1==q)) update_jeu(Inter[q],process);
                      } else if (process.parallelisation->is_master_cpu()) std::cout << "La fonction donn�e n'est pas sous la bonne forme : 0;0;0" << std::endl;
                    } else if (Inter[q].comp=="Contact" or Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique"){
                      if (process.parallelisation->is_master_cpu()) std::cout << "Interface de type " << Inter[q].comp << " donner le nouveau coefficient de frottement " << Inter[q].param_comp->coeffrottement << std::endl;
                      if (process.parallelisation->is_master_cpu()) {
                        if (f) {
                          getline(f,str);
                          std::cout << str << std::endl;
                        } else {
                          getline(std::cin,str);
                          process.affichage->command_file="";
                        }
                      }
                      if (process.parallelisation->is_multi_cpu()) MPI_Bcast(str,0);
                      Inter[q].param_comp->coeffrottement=atof(str.c_str());
                    } else {
                      if (process.parallelisation->is_master_cpu()) std::cout << "Interface de type " << Inter[q].comp << " non modifiable" << std::endl;
                    }
                  } else {
                    if (process.parallelisation->is_master_cpu()) std::cout << "mauvais numero d interface" << std::endl;
                  }
                  } else {
                    if (process.parallelisation->is_master_cpu()) std::cout << "donner un numero d interface" << std::endl;
                  }
                }
                if(v[0]=="calcul") {
                    process.affichage->fichiers_paraview_inter_crees=0;
                    process.affichage->fichiers_paraview_sst_crees=0;
                    process.reprise_calcul=2;

//                     if(process.nom_calcul=="incr") multiscale_iterate_incr(S,SubS, Inter,SubI, process, Global, CL);
//                     else if(process.nom_calcul=="latin") multiscale_iterate_latin(S,SubS, Inter,SubI, process, Global, CL);


                }
            }
        }
    }
}


