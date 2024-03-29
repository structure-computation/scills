#include "assignation_material_properties_Interface.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"
#include "../DEFINITIONS/ParallelisationData.h"


void assignation_materials_property_INTER(DataUser &data_user, Vec<Interface> &Inter, Process &process, FieldStructureUser &field_structure_user){
    process.print("\t Assignation de materiau particulier (Ex : contact) aux Interfaces");
    Vec<InterCarac > propinter;
    BasicVec<BasicVec<TYPEREEL > > link_prop_temp;
    read_propinter(propinter,data_user, link_prop_temp);

    modif_inter(Inter,propinter,process,data_user);
}


void read_propinter(Vec<InterCarac> &propinter,const DataUser &data_user, BasicVec<BasicVec<TYPEREEL> > &link_prop_temp) {
    unsigned nbliaisons = data_user.links_vec.size();
    //propinter.resize(nbliaisons);   // Fais dans le Process::read_data_user
    for(unsigned i=0;i<nbliaisons;++i) {
        const DataUser::Json_links &link = data_user.links_vec[i];
        propinter[i].id = link.id_in_calcul;
        switch(link.type_num){
            case 0:
                propinter[i].type = "parfait";
                propinter[i].comp = "Parfait";
                break;
            case 1:
                propinter[i].type = "elastique";
                propinter[i].comp = "Elastique";
                break;
            case 2:
                propinter[i].type = "contact parfait";
                propinter[i].comp="Contact Parfait";
                break;
            case 3:
                propinter[i].type = "parfait cassable";
                propinter[i].comp = "Parfait Cassable";
                break;
            case 4:
                propinter[i].type = "elastique cassable";
                propinter[i].comp = "Elastique Cassable";
                break;
            case 5:
                propinter[i].type = "cohesive";
                propinter[i].comp = "Cohesive";
                break;
            default:
                std::cerr  << "comportement d'interface non reconnu : type_num = " << link.type_num << std::endl;
                assert(0);
        }
        /// Caracteristique du prechargement
        propinter[i].Ep_Type        = link.Ep_type;     /// Type de condition interne (epaisseur/precharge impose(e)/normal(e))
        propinter[i].Ep_n           = link.Ep_n;        /// jeux ou epaisseur normal
        propinter[i].Ep_x           = link.Ep_x;        /// jeux ou epaisseur selon x
        propinter[i].Ep_y           = link.Ep_y;        /// jeux ou epaisseur selon y
        propinter[i].Ep_z           = link.Ep_z;        /// jeux ou epaisseur selon z
        propinter[i].Preload_n      = link.Preload_n;   /// precharge normale
        propinter[i].Preload_x      = link.Preload_x;   /// precharge selon x
        propinter[i].Preload_y      = link.Preload_y;   /// precharge selon y
        propinter[i].Preload_z      = link.Preload_z;   /// precharge selon z
        
        /// Caracteristiques pour le contact avec frottement
        propinter[i].f              = link.f;           /// coefficient de frottement
        
        /// Caracteristiques pour le comportement cassable
        propinter[i].Fcr_n          = link.Fcr_n;       /// limite en rupture normale
        propinter[i].Fcr_t          = link.Fcr_t;       /// limite en rupture tangentielle
        
        /// Caracteristiques pour le comportement elastique
        propinter[i].Kn             = link.Kn;          /// Raideur normale (globale ou en traction) de l'interface
        propinter[i].Kt             = link.Kt;          /// Raideur tangentielle de l'interface
        propinter[i].Knc            = link.Knc;         /// Raideur normale en compression de l'interface
        
        /// Caracteristiques pour le comportement plastique
        propinter[i].Rop            = link.Rop;         /// Contrainte initiale de plastification
        propinter[i].kp             = link.kp;          /// Multiplicateur de la loi d'ecrouissage
        propinter[i].mp             = link.np;          /// Exposant de la loi d'ecrouissage
        
        /// Caracteristiques pour les interfaces cohesives (type mesomodele)
        propinter[i].Yo             = link.Yo;          /// Valeur minimale d'effort thermodynamique equivalent pour amorcer l'endommagement
        propinter[i].Yc             = link.Yc;          /// Valeur critique d'effort thermodynamique equivalent pour l'endommagement
        propinter[i].alpha          = link.alpha;       /// exposant dans le calcul de l'effort thermodynamique equivalent
        propinter[i].gamma          = link.gamma;       /// couplage normal/tangentiel dans le calcul de l'effort thermodynamique equivalent
        propinter[i].n              = link.n;           /// coefficient dans le calcul de l'endommagement associe a un effort thermodynamique
    }
    /*
    unsigned nbliaisons = data_user.behaviour_links.size();
    propinter.resize(nbliaisons);
    link_prop_temp.resize(nbliaisons);
    for(unsigned i=0;i<nbliaisons;++i) {
        propinter[i].id = data_user.behaviour_links[i].id;
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
            std::cerr  << "comportement d'interface non reconnu" << std::endl;
            assert(0);
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
    }*/
}


void modif_inter(Vec<Interface> &Inter, Vec<InterCarac> &propinter, Process &process,const DataUser &data_user) {
    //assignation des proprietes materiau des interfaces
    for(unsigned i=0;i<propinter.size();i++) {
        LMT::Vec<unsigned> inter_select;
        //reperage par le champ id_link de l'Interface
        for(unsigned q=0;q<Inter.size();q++){
            if (Inter[q].type=="Int" and Inter[q].id_link == propinter[i].id) {
                inter_select.push_back(q);
            }
        }
        
        //assignation des proprietes aux interfaces selectionnees
        for(unsigned j=0;j<inter_select.size();j++) {
            unsigned q=inter_select[j];
            //if (process.parallelisation->is_master_cpu()) std::cout << "\t  Interface modifiee : " << Inter[q].num << std::endl;
            Inter[q].comp=propinter[i].comp;
            if (propinter[i].type=="contact_sst" or propinter[i].type=="contact_box") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
            }
            if (propinter[i].type=="contact_jeu_sst" or propinter[i].type=="contact_jeu_box" or propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" or propinter[i].type=="contact_ep") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
                //if (propinter[i].type=="jeu_impose_sst" or propinter[i].type=="jeu_impose_box" ) Inter[q].param_comp->nbpastempsimpos=propinter[i].nbpastempsimpos;
                /*std::vector<Ex> symbols;
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
                *//*
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
                }*/
            }
            if (propinter[i].type=="contact_jeu_physique") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
                //le champ jeu est assigne quand on fait la correspondances des elements des deux maillages d interfaces, comme on a la distance g1 g2 dans op_inter.h
            }
            if (propinter[i].type=="discrete") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
                //Inter[q].param_comp->Gcrit=propinter[i].Gcrit;
            }
            if (propinter[i].type=="cohesive") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
                //Inter[q].param_comp->param_damage=propinter[i].param_damage;
            }
            if (propinter[i].type=="breakable") {
                //Inter[q].param_comp->coeffrottement=propinter[i].coeffrottement;
                //Inter[q].param_comp->jeu.set(0.);
                //Inter[q].param_comp->Gcrit=propinter[i].Gcrit;
            }
            if (propinter[i].type=="contact_ep" or propinter[i].type=="parfait") {
                //Inter[q].param_comp->coeffrottement=0.;
                //Inter[q].param_comp->f_coeffrottement.set(0.);
                //Inter[q].param_comp->jeu.set(0.);
                /*std::vector<Ex> symbols;
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
                }*/
            }
        }
    }
    /*
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
            //if (process.parallelisation->is_master_cpu()) std::cout << "\t  Interface modifiee : " << Inter[q].num << std::endl;
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
                //le champ jeu est assigne quand on fait la correspondances des elements des deux maillages d interfaces, comme on a la distance g1 g2 dans op_inter.h
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
    }//*/
}