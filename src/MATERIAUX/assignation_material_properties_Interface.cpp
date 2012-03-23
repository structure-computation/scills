#include "assignation_material_properties_Interface.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"


void assignation_materials_property_INTER(DataUser &data_user, Vec<Interface> &Inter, Param &process, FieldStructureUser &field_structure_user){
    if (process.rank == 0) std::cout << "\t Assignation de materiau particulier (Ex : contact) aux Interfaces" << std::endl;
    Vec<InterCarac > propinter;
    BasicVec<BasicVec<TYPEREEL > > link_prop_temp;
    read_propinter(propinter,data_user, link_prop_temp);
    PRINT("lecture propriete interface ok");

    modif_inter(Inter,propinter,process,data_user);
    //idem pour les group_interfaces (pour GPU)
    field_structure_user.assign_link_id_to_group_interfaces(data_user);
    field_structure_user.assign_link_properties_to_group_interfaces(data_user,link_prop_temp);
}


void read_propinter(Vec<InterCarac> &propinter,const DataUser &data_user, BasicVec<BasicVec<TYPEREEL> > &link_prop_temp) {
    unsigned nbliaisons = data_user.behaviour_links.size();
    propinter.resize(nbliaisons);
    link_prop_temp.resize(nbliaisons);
    for(unsigned i=0;i<nbliaisons;++i) {
        propinter[i].id = data_user.behaviour_links[i].id;
        //         PRINT(data_user.behaviour_links[i].type_num);
        if(data_user.behaviour_links[i].type_num == 0){   //parfaite
            propinter[i].type = "parfait";
            propinter[i].comp="Parfait";
        }else if(data_user.behaviour_links[i].type_num == 1){   //elastique
            propinter[i].type = "elastique";
            propinter[i].comp="Parfait";
        }else if(data_user.behaviour_links[i].type_num == 2){
            //             propinter[i].type = "contact_jeu_sst";
            //             propinter[i].comp="Contact_jeu"; 
            propinter[i].type = "contact_ep";
            propinter[i].comp="Contact_ep";
        }else if(data_user.behaviour_links[i].type_num == 3){
            propinter[i].type = "cohesive";
            propinter[i].comp="Cohesive";
        }else if(data_user.behaviour_links[i].type_num == 4){
            propinter[i].type = "breakable";
            propinter[i].comp="Breakable";
        }else{
            std::cout  << "comportement d'interface non reconnu" << std::endl;  assert(0);
        }
        
        
        std::vector<Ex> symbols;
        if (DIM==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (DIM==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }
        
        if(data_user.options.Multiresolution_on==1){
            //ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
                symbols.push_back(var.c_str());
            }
        }
        
        link_prop_temp[i].resize(data_user.behaviour_links[i].link_prop.size());
        for(int i_prop=0; i_prop<data_user.behaviour_links[i].link_prop.size(); i_prop++){
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_links[i].link_prop[i_prop],symbols);
            Ex::MapExNum var_temp;
            for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var_temp[symbols[d2]]= 0.;
            }
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            link_prop_temp[i][i_prop] = (TYPEREEL) expr_temp.subs_numerical(var_temp);           
        }
        propinter[i].f_coeffrottement = data_user.behaviour_links[i].link_prop[0];      /// coeff frottement analytique
        propinter[i].coeffrottement = link_prop_temp[i][0];                             /// coeff frottement
        propinter[i].jeu = data_user.behaviour_links[i].link_prop[1];                   /// jeux ou epaisseur negative        
        propinter[i].Gcrit = link_prop_temp[i][7];                                      /// limite en rupture    
        propinter[i].f_R = data_user.behaviour_links[i].link_prop[3];                   /// coeff frottement analytique
        //PRINT(link_prop_temp[i]);
    }
}


void modif_inter(Vec<Interface> &Inter, Vec<InterCarac> &propinter, Param &process,const DataUser &data_user) {
    //allocation de la memoire pour les parametres de comportement d'interface
    for(unsigned q=0;q<Inter.size();++q) {
        #ifdef PRINT_ALLOC
        total_allocated[ typeid(PARAM_COMP_INTER).name() ] += sizeof(PARAM_COMP_INTER);
        #endif
        Inter[q].param_comp = new PARAM_COMP_INTER(Inter[q].side[0].nodeeq.size());
        Inter[q].param_comp->jeu.resize(Inter[q].side[0].nodeeq.size()*DIM,0.);
        //  std::cout << "DEBUG : " << Inter[q].side[0].nodeeq.size()*DIM << std::endl;
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
        
        //assignation des proprietes aux interfaces selectionnees
        for(unsigned j=0;j<inter_select.size();j++) {
            unsigned q=inter_select[j];
            //if (process.rank == 0) std::cout << "\t  Interface modifiee : " << Inter[q].num << std::endl;
            Inter[q].comp=propinter[i].comp;
            if (propinter[i].type=="contact_sst" or propinter[i].type=="contact_box") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
            }
            if (propinter[i].type=="contact_jeu_sst" or propinter[i].type=="contact_jeu_box" or propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" or propinter[i].type=="contact_ep") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                if (propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" ) Inter[q].param_comp->nbpastempsimpos=propinter[i].nbpastempsimpos;
                std::vector<Ex> symbols;
                if (DIM==2) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                } else if (DIM==3) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                    symbols.push_back("z");
                }
                
                //ajout des variables de multiresolution aux symboles
                if(data_user.options.Multiresolution_on==1){
                    for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                        Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
                        symbols.push_back(var.c_str());
                    }
                }
                
                Vec<Sc2String> jeu_cut=tokenize(propinter[i].jeu,';');
                Inter[q].param_comp->fcts_spatiales=propinter[i].jeu;
                if (jeu_cut.size() == 1) {//Jeu normal
                    Ex expr;
                    expr = read_ex(propinter[i].jeu.c_str(),symbols);
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        TYPEREEL data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        if(data_user.options.Multiresolution_on==1){
                            //evaluation des variables de multiresolution
                            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                        }
                        
                        data = (TYPEREEL)expr.subs_numerical(var);
                        Inter[q].param_comp->jeu[range(DIM*k,DIM*(k+1))]=data*Inter[q].side[0].neq[range(DIM*k,DIM*(k+1))];
                    }
                } else {//Jeu complet
                    Vec<Ex> expr;
                    expr.resize(DIM);
                    
                    for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                        expr[d2] = read_ex(jeu_cut[d2],symbols);
                    }
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        Vec<TYPEREEL,DIM> data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        }
                        if(data_user.options.Multiresolution_on==1){
                            //evaluation des variables de multiresolution
                            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                        }
                        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            data[d2] = (TYPEREEL)expr[d2].subs_numerical(var);
                        Inter[q].param_comp->jeu[range(DIM*k,DIM*(k+1))]=data;
                        
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
            if (propinter[i].type=="breakable") {
                Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                Inter[q].param_comp->jeu.set(0.);
                Inter[q].param_comp->Gcrit=propinter[i].Gcrit;
            }
            if (propinter[i].type=="contact_ep" or propinter[i].type=="parfait") {
                Inter[q].param_comp->coeffrottement=0.;
                Inter[q].param_comp->f_coeffrottement.set(0.);
                Inter[q].param_comp->jeu.set(0.);
                std::vector<Ex> symbols;
                if (DIM==2) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                } else if (DIM==3) {
                    symbols.push_back("x");
                    symbols.push_back("y");
                    symbols.push_back("z");
                }
                //ajout des variables de multiresolution aux symboles
                if(data_user.options.Multiresolution_on==1){
                    for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                        Sc2String var="V"; var<<i_par; //nom de la variable de multiresolution
                        symbols.push_back(var.c_str());
                    }
                }
                Vec<Sc2String> jeu_cut=tokenize(propinter[i].jeu,';');
                Vec<Sc2String> f_coeffrottement_cut=tokenize(propinter[i].f_coeffrottement,';');
                
                // coefficient de frottement
                Inter[q].param_comp->fcts_spatiales=propinter[i].f_coeffrottement;
                TYPEREEL sum_data=0;
                if (f_coeffrottement_cut.size() == 1) {//Jeu normal
                    Ex expr;
                    expr = read_ex(propinter[i].f_coeffrottement.c_str(),symbols);
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        TYPEREEL data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        if(data_user.options.Multiresolution_on==1){
                            //evaluation des variables de multiresolution
                            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                        }
                        data = (TYPEREEL)expr.subs_numerical(var);
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
                        TYPEREEL data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        if(data_user.options.Multiresolution_on==1){
                            //evaluation des variables de multiresolution
                            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                        }
                        data = (TYPEREEL)expr.subs_numerical(var);
                        Inter[q].param_comp->jeu[range(DIM*k,DIM*(k+1))]=data*Inter[q].side[0].neq[range(DIM*k,DIM*(k+1))];
                    }
                } else {//Jeu complet
                    Vec<Ex> expr;
                    expr.resize(DIM);
                    
                    for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                        expr[d2] = read_ex(jeu_cut[d2],symbols);
                    }
                    for(unsigned k=0;k<Inter[q].side[0].nodeeq.size();++k) {
                        Vec<TYPEREEL,DIM> data;
                        Ex::MapExNum var;
                        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                            var[symbols[d2]]= Inter[q].side[0].nodeeq[k][d2];
                        }
                        if(data_user.options.Multiresolution_on==1){
                            //evaluation des variables de multiresolution
                            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
                        }
                        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
                            data[d2] = (TYPEREEL)expr[d2].subs_numerical(var);
                        Inter[q].param_comp->jeu[range(DIM*k,DIM*(k+1))]=data;
                    }
                }
            }
        }
    }
}