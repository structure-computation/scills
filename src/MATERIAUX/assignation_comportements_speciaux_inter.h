using namespace LMT;
using namespace std;
#include "definition_PARAM_MPI.h"

/** \ingroup Materiaux
\brief Assignation des propriétés aux interfaces
 
On alloue tout d'abord la mémoire pour les paramètres matériaux (cf. PARAM_COMP_INTER) de toutes les interfaces.
 
On boucle ensuite sur les propriétés d'interfaces (Matprop) et selon le type (box ou sst) on sélectionne les interfaces correspondantes.
 
Enfin pour ces interfaces sélectionnées seulement, on modifie le champ Interface::comp, on assigne ensuite les paramètres aux interfaces.
 
Pour le jeu, en utilisant codegen et les noeuds équivalents de l'interface, on crée le champ jeu qui sera par la suite utilisé dans l'étape locale. Si l'utilisateur renseigne une seule valeur, on considère que c'est un jeu selon la normale, sinon il faut renseigner les valeurs selon x, y (et z) du jeu. 
*/
template<class TV2, class TV4, class TV1>
void modif_inter(TV2 &Inter, TV4 &propinter, TV1 &S,Param &process) {
    //allocation de la memoire pour les parametres de comportement d'interface
    for(unsigned q=0;q<Inter.size();++q) {
#ifdef PRINT_ALLOC
        total_allocated[ typeid(PARAM_COMP_INTER).name() ] += sizeof(PARAM_COMP_INTER);
#endif
        Inter[q].param_comp = new PARAM_COMP_INTER(Inter[q].side[0].nodeeq.size());
        Inter[q].param_comp->jeu.resize(Inter[q].side[0].nodeeq.size()*TV2::template SubType<0>::T::dim);
    }

    //assignation des proprietes materiau des interfaces
    for(unsigned i=0;i<propinter.size();i++) {
        Vec<unsigned> inter_select;
        //reperage par le champ id_link de l'Interface
        for(unsigned q=0;q<Inter.size();q++){
            if (Inter[q].type=="Int" and Inter[q].id_link == propinter[i].id) {
                inter_select.push_back(q);
            }
        }
        
//         if (propinter[i].type=="contact_sst" or propinter[i].type=="contact_jeu_sst" or propinter[i].type=="contact_jeu_physique" or propinter[i].type=="jeu_impose_sst") {
//             //reperage par le numero des ssts voisines :
//             Vec<unsigned,2> numsst = propinter[i].num_sst;
//             Vec<Vec<unsigned>,2> numinter;
//             for(unsigned k=0;k<2;k++)
//                 for(unsigned j=0;j<S[numsst[k]].edge.size();j++)
//                     numinter[k].push_back(S[numsst[k]].edge[j].internum);
// 
//             for(unsigned k=0;k<numinter[0].size();k++)
//                 if(find(numinter[1],LMT::_1==numinter[0][k])==1) {
//                 if(find(process.multi_mpi->listinter,LMT::_1==numinter[0][k])){
//                     inter_select.push_back(numinter[0][k]);
//                     break;
//                 }
//                 }
//         } else {
//             for(unsigned q=0;q<Inter.size();q++)
//                 if (Inter[q].type=="Int") {
//                 if (find(process.multi_mpi->listinter,LMT::_1==q)) {
//                     if (pt_in_box(Inter[q].G,propinter[i].box)==1)
//                         inter_select.push_back(q);
//                 }
//                 }
//         }

        //assignation des proprietes aux interfaces selectionnees
        for(unsigned j=0;j<inter_select.size();j++) {
            unsigned q=inter_select[j];
            if (process.rank == 0) std::cout << "\t  Interface modifiee : " << Inter[q].num << std::endl;
            Inter[q].comp=propinter[i].comp;
            if (propinter[i].type=="contact_sst" or propinter[i].type=="contact_box") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
            }
            if (propinter[i].type=="contact_jeu_sst" or propinter[i].type=="contact_jeu_box" or propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" or propinter[i].type=="contact_ep") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                if (propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" ) Inter[q].param_comp->nbpastempsimpos=propinter[i].nbpastempsimpos;
                typedef typename TV2::template SubType<0>::T::T T;
                std::vector<Ex> symbols;
                if (TV2::template SubType<0>::T::dim==2) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                } else if (TV2::template SubType<0>::T::dim==3) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                    symbols.push_back("z");
                }
                Vec<string> jeu_cut=tokenize(propinter[i].jeu,';');
                Inter[q].param_comp->fcts_spatiales=propinter[i].jeu;
                if (jeu_cut.size() == 1) {//Jeu normal
                    Ex expr;
                    expr = read_ex(propinter[i].jeu.c_str(),symbols);
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        T data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        data = (T)expr.subs_numerical(var);
                        Inter[q].param_comp->jeu[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))]=data*Inter[q].side[0].neq[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))];
                    }
                } else {//Jeu complet
                    Vec<Ex> expr;
                    expr.resize(TV2::template SubType<0>::T::dim);

                    for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                        expr[d2] = read_ex(jeu_cut[d2],symbols);
                    }
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        Vec<T,TV2::template SubType<0>::T::dim> data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        }
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            data[d2] = (T)expr[d2].subs_numerical(var);
                        Inter[q].param_comp->jeu[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))]=data;

                    }
                }
            }
            if (propinter[i].type=="contact_jeu_physique") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                //le champ jeu est assigné quand on fait la correspondances des éléments des deux maillages d interfaces, comme on a la distance g1 g2 dans op_inter.h
            }
            if (propinter[i].type=="discrete") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                Inter[q].param_comp->Gcrit=propinter[i].Gcrit;
            }
            if (propinter[i].type=="cohesive") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                Inter[q].param_comp->param_damage=propinter[i].param_damage;
            }
            if (propinter[i].type=="contact_ep" or propinter[i].type=="parfait") {
                std::cout << "propinter[i].type = " << propinter[i].type << std::endl;
                
                std::cout << "  propinter[i].f_coeffrottement = " << propinter[i].f_coeffrottement << std::endl;
                std::cout << "  propinter[i].jeu = " << propinter[i].jeu << std::endl;
                
                Inter[q].param_comp->coeffrottement=0.;
                Inter[q].param_comp->f_coeffrottement.set(0.);
                Inter[q].param_comp->jeu.set(0.);
                typedef typename TV2::template SubType<0>::T::T T;
                std::vector<Ex> symbols;
                if (TV2::template SubType<0>::T::dim==2) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                } else if (TV2::template SubType<0>::T::dim==3) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                    symbols.push_back("z");
                }
                Vec<string> jeu_cut=tokenize(propinter[i].jeu,';');
                Vec<string> f_coeffrottement_cut=tokenize(propinter[i].f_coeffrottement,';');

                // coefficient de frottement
                Inter[q].param_comp->fcts_spatiales=propinter[i].f_coeffrottement;
                T sum_data=0;
                if (f_coeffrottement_cut.size() == 1) {//Jeu normal
                    Ex expr;
                    expr = read_ex(propinter[i].f_coeffrottement.c_str(),symbols);
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        T data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        data = (T)expr.subs_numerical(var);
                        Inter[q].param_comp->f_coeffrottement[k]=data;
                        sum_data += data;
                    }
                    Inter[q].param_comp->coeffrottement=sum_data/Inter[q].side[0].nodeeq.size();
                }

                // defaut de forme
                Inter[q].param_comp->fcts_spatiales=propinter[i].jeu;
                if (jeu_cut.size() == 1) {//Jeu normal
                    Ex expr;
                    expr = read_ex(propinter[i].jeu.c_str(),symbols);
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        T data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        data = (T)expr.subs_numerical(var);
                        Inter[q].param_comp->jeu[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))]=data*Inter[q].side[0].neq[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))];
                    }
                } else {//Jeu complet
                    Vec<Ex> expr;
                    expr.resize(TV2::template SubType<0>::T::dim);

                    for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                        expr[d2] = read_ex(jeu_cut[d2],symbols);
                    }
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        Vec<T,TV2::template SubType<0>::T::dim> data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        }
                        for(unsigned d2=0;d2<TV2::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            data[d2] = (T)expr[d2].subs_numerical(var);
                        Inter[q].param_comp->jeu[range(TV2::template SubType<0>::T::dim*k,TV2::template SubType<0>::T::dim*(k+1))]=data;

                    }
                }
            }

        }

    }
}
