#include "Interface.h"
#include "Boundary.h"

//---------------------------- Attributs statiques ----------------------------//
Sc2String Interface::type_ext = "Ext";                                  /// Nom pour le type interface exterieure (CL)
Sc2String Interface::type_int = "Int";                                  /// Nom pour le type interface interieure (liaison)
Sc2String Interface::comp_parfait = "Parfait";                          /// Nom pour le comportement parfait
Sc2String Interface::comp_contact_parfait = "Contact";                  /// Nom pour le comportement contact de type parfait
Sc2String Interface::comp_contact_elastique =  "Contact elastique";     /// Nom pour le comportement contact de type elastique
Sc2String Interface::comp_elastique = "Elastique";                      /// Nom pour le comportement elastique
Sc2String Interface::comp_cohesive = "Cohesive";                        /// Nom pour le comportement cohesif
Sc2String Interface::comp_cassable_parfait = "Parfait Cassable";        /// Nom pour le comportement parfait cassable
Sc2String Interface::comp_cassable_elastique = "Elastique Cassable";    /// Nom pour le comportement elastique cassable
Sc2String Interface::comp_deplacement = "depl";                         /// Nom pour le comportement deplacement impose
Sc2String Interface::comp_deplacement_normal = "depl_normal";           /// Nom pour le comportement deplacement normal
Sc2String Interface::comp_deplacement_nul = "depl_nul";                 /// Nom pour le comportement deplacement nul
Sc2String Interface::comp_vitesse = "vit";                              /// Nom pour le comportement vitesse imposee
Sc2String Interface::comp_vitesse_normale = "vit_normale";              /// Nom pour le comportement vitesse normale
Sc2String Interface::comp_vitesse_nulle = "vit_nulle";                  /// Nom pour le comportement vitesse nulle
Sc2String Interface::comp_effort = "effort";                            /// Nom pour le comportement effort impose
Sc2String Interface::comp_effort_normal = "effort_normal";              /// Nom pour le comportement effort normal (pression)
Sc2String Interface::comp_symetrie = "sym";                             /// Nom pour le comportement symetrie
Sc2String Interface::comp_periodique = "periodique";                    /// Nom pour le comportement periodique
//-----------------------------------------------------------------------------//

Interface::Side::Side(){
    t.resize(1);
    mesh=NULL;
}

/// fonction projecteur macro . La définition d'un projecteur macro et micro sous forme d'une fonction permet de ne pas construire et stocker de nouvel opérateur et prend autant de temps que l'utilisation d'une matrice.
Vector Interface::Side::PM(Vector &f) {
    Vector fM;
    fM.resize(f.size());
    fM.set(0.);
    for(unsigned i=0;i<eM.nb_cols();i++){
        fM+=dot(eM.col(i),f)*eM.col(i);
    }
    return fM;
}
/// fonction projecteur micro
Vector Interface::Side::Pm(Vector &f) {
    Vector fm;
    fm.resize(f.size());
    fm.set(0.);
    for(unsigned i=0;i<eM.nb_cols();i++){
        fm+=dot(eM.col(i),f)*eM.col(i);
    }
    return f-fm;
}
/// fonction projecteur normal
Vector Interface::Side::Pn(Vector &f) {
    Vector fn;
    fn.resize(f.size());
    fn.set(0.);
    for(unsigned i=0;i<nodeeq.size();++i){
        LMT::Vec<int,DIM> rep=range(i*DIM,(i+1)*DIM);
        fn[rep]=dot(neq[rep],f[rep])*neq[rep];
    }
    return fn;
}
/// fonction projecteur tangentiel
Vector Interface::Side::Pt(Vector &f) {
    return f-Pn(f);
}

void Interface::Side::Time::allocations(unsigned sizenodeeq,bool endommageable){
    F.resize(sizenodeeq);
    F.set(0.0);
    Wp.resize(sizenodeeq);
    Wp.set(0.0);
    W.resize(sizenodeeq);
    W.set(0.0);
    oldW.resize(sizenodeeq);
    oldW.set(0.0);
    Fchap.resize(sizenodeeq);
    Fchap.set(0.0);
    Wchap.resize(sizenodeeq);
    Wchap.set(0.0);
    Wpchap.resize(sizenodeeq);
    Wpchap.set(0.0);
    WtildeM.resize(sizenodeeq);
    WtildeM.set(0.0);
    oldF.resize(sizenodeeq);
    oldF.set(0.0);
    oldWp.resize(sizenodeeq);
    oldWp.set(0.0);
    
    if(endommageable){
        d.resize(sizenodeeq/DIM);
        d.set(0.0);
    }
}

void Interface::allocate(unsigned int nbpastemps){
    for(unsigned j=0;j<side.size();++j){
        side[j].t.resize(nbpastemps+1);
        //if(process.recopie_t_post==1){
        //    side[j].t=side[j].t_post;
        //}
        //else{
            const unsigned sizenodeeq = side[j].nodeeq.size()*DIM;
            for(unsigned pt=0;pt<(nbpastemps+1);pt++){
                /// Allocation des vecteurs de prechargement
                Ep_imposee.resize(sizenodeeq);
                Ep_imposee.set(0.0);
                old_Ep_imposee.resize(sizenodeeq);
                old_Ep_imposee.set(0.0);
                precharge.resize(sizenodeeq);
                precharge.set(0.0);
                old_precharge.resize(sizenodeeq);
                old_precharge.set(0.0);
                /// Allocation des vecteurs d'elasticite
                if(comp == Interface::comp_elastique or
                    comp == Interface::comp_cassable_elastique or
                    comp == Interface::comp_contact_elastique or
                    comp == Interface::comp_cohesive){
                    Ep_elastique.resize(sizenodeeq);
                    Ep_elastique.set(0.0);
                    old_Ep_elastique.resize(sizenodeeq);
                    old_Ep_elastique.set(0.0);
                }
                /// Si le comportement peut etre modifie durant le calcul
                if(comp == Interface::comp_cassable_parfait or 
                    comp == Interface::comp_cassable_elastique or 
                    comp == Interface::comp_cohesive){
                    comportement.resize(sizenodeeq);
                    comportement.set(false);
                }
                side[j].t[pt].allocations(sizenodeeq,(matprop == 0)?false:matprop->degradable);
            }
        //}
    }
}


///suppression des interfaces
void Interface::free(){
    vois.free();
    repddl.free();

    #ifdef PRINT_ALLOC
    if (side[0].mesh != NULL) total_allocated[ typeid(typename Interface::TMESH).name() ] -= sizeof(typename Interface::TMESH);
    #endif

    /// Les maillages sont les memes des 2 cotes...
    if (side.size() != 0){
        if (side[0].mesh != NULL) delete side[0].mesh;
        side[0].mesh=NULL;
        side.free();
    }
    /// matprop est gere par le vecteur de SstCarac
}

/// Initialisation des proprietes de l'interface
void Interface::init(unsigned pt){
    const LMT::Vec<Point> &nodes = side[0].nodeeq;
    const int nb_nodes = nodes.size();
    Ep_imposee.resize(nb_nodes*DIM,0.0);
    if(old_Ep_imposee.size()!=Ep_imposee.size()){
        old_Ep_imposee.resize(nb_nodes*DIM,0.0);
    }
    Vector preload_temp;
    preload_temp.resize(nb_nodes*DIM,0.0);
    coeffrottement_vec.resize(nb_nodes,0.0);
    coeffrottement = 0;
    if(matprop == 0){
        return;
    }
    Ex::MapExNum values = InterCarac::inter_materials_parameters.getParentsValues();             /// Recuperation des parametres temporels et de multiresolution
    for(unsigned i=0;i<nb_nodes;++i){
        for(unsigned i_dir=0;i_dir<DIM;++i_dir){
            values[Boundary::CL_parameters.main_parameters[i_dir]->self_ex] = nodes[i][i_dir];  /// Chargement des coordonnees du point (main_parameters pointe vers x, y et z)
        }
        LMT::Vec<int,DIM> rep = range(i*DIM,(i+1)*DIM);
        //PRINT(matprop->Ep_Type);
        switch(matprop->Ep_Type){
            case 0:
                /// Jeu normal
                preload_temp[rep] = matprop->Ep_n.updateValue(values)*side[0].neq[rep];
                break;
            case 1:
                /// Jeu impose
                preload_temp[rep[0]] = matprop->Ep_x.updateValue(values);
                preload_temp[rep[1]] = matprop->Ep_y.updateValue(values);
                #if DIM == 3
                preload_temp[rep[2]] = matprop->Ep_z.updateValue(values);
                #endif
                break;
            case 2:
                /// Precharge normal
                preload_temp[rep] = matprop->Preload_n.updateValue(values)*side[0].neq[rep];
                break;
            case 3:
                /// Precharge imposee
                preload_temp[rep[0]] = matprop->Preload_x.updateValue(values);
                preload_temp[rep[1]] = matprop->Preload_y.updateValue(values);
                #if DIM == 3
                preload_temp[rep[2]] = matprop->Preload_z.updateValue(values);
                #endif
                break;
            default:
                break;
        }
        coeffrottement_vec[i] = matprop->f.updateValue(values);     /// Evaluation des composantes du frottement
        coeffrottement += coeffrottement_vec[i];                    /// incrementation du coefficient de frottement global
    }
    switch(matprop->Ep_Type){
        case 0:
            /// Jeu normal
            Ep_imposee = side[0].Pn(preload_temp);
            break;
        case 1:
            /// Jeu impose
            Ep_imposee = preload_temp;
            break;
        case 2:
            /// Precharge normal
            precharge = side[0].Pn(preload_temp)/measure;
//             PRINT(preload_temp/measure);
//             PRINT(side[0].Pn(preload_temp));
//             PRINT(side[0].Pn(preload_temp)/measure);
//             PRINT(measure);
            break;
        case 3:
            /// Precharge imposee
            precharge = preload_temp/measure;
            break;
        default:
            break;
    }
    coeffrottement /= nb_nodes;
    if(matprop->Ep_Type == 0 or matprop->Ep_Type == 1){
        PRINT(Ep_imposee[LMT::range(0,DIM)]);
    }
    else if(matprop->Ep_Type == 2 or matprop->Ep_Type == 3){
        PRINT(precharge[LMT::range(0,DIM)]);
    }
}

void Interface::affiche(){
    std::cout << "------------------- interface ------------------------" << std::endl;
    PRINT(id);
    PRINT(type);
    PRINT(comp);
    PRINT(matprop);
    PRINT(refCL);
    PRINT(edge_id);
    PRINT(id_bc);
    PRINT(id_link);
    PRINT(vois);
    std::cout << "------------------- fin interface ------------------------" << std::endl;
}

void Interface::read_data_user(int index, const DataUser& data_user, const GeometryUser &geometry_user)
{/*
    const DataUser::Json_interfaces &interface = data_user.interfaces_vec[index];
    id = interface.id_in_calcul;
    num = interface.id_in_calcul;
    id_link = data_user.find_links_index(interface.link_id);
    
    std::istringstream s(interface.adj_num_group);
    int id_sst_0,id_sst_1;
    s>>id_sst_0;
    s>>id_sst_1;
    id_sst_0 = data_user.find_pieces_index(id_sst_0);
    id_sst_1 = data_user.find_pieces_index(id_sst_1);
*/}

int Interface::get_type_elem() const {
    int type_elem=0;
    if (type == type_ext and (comp == comp_deplacement or comp == comp_vitesse or comp == comp_deplacement_nul or comp == comp_vitesse_nulle)){
        type_elem = 0;
    } else if (type == type_ext and (comp == comp_effort)){
        type_elem = 1;
    } else if (type == type_ext and (comp == comp_symetrie)){
        type_elem = 2;
    } else if (type == type_ext and (comp == comp_deplacement_normal or comp == comp_vitesse_normale)){
        type_elem = 3;
    } else if (type == type_int and (comp == comp_parfait)){
        type_elem = 4;
    } else if (type == type_int and (comp=="Contact_jeu_physique" or comp == comp_contact_parfait or comp == comp_contact_elastique) ){
        type_elem = 5;
    } else if (type == type_int and (comp=="Jeu_impose")){
        type_elem = 6;
    } else if (type == type_ext and (comp == comp_periodique)){
        type_elem = 7;
    } else {
        type_elem = 8;
    }
    return type_elem;
}

Interface::Interface() {matprop = 0;id_link = -1;}
Interface::~Interface() {free();}



/**********************************************************************************************
*                                                                                             *
*                                         NodalState                                          *
*                                                                                             *
**********************************************************************************************/

Interface::NodalState::NodalState(Interface &I,unsigned pt,Scalar dt_):
interface(I),
matprop(I.matprop),
i_time(pt),
dt(dt_)
{
    Side::Time& side1 = interface.side[0].t[i_time];
    Side::Time& side2 = interface.side[1].t[i_time];
    Side::Time& old_side1 = interface.side[0].t[i_time-1];
    Side::Time& old_side2 = interface.side[1].t[i_time-1];
    LMT::Vec<unsigned> &list1 = interface.side[0].ddlcorresp;
    LMT::Vec<unsigned> &list2 = interface.side[1].ddlcorresp;
    const unsigned nb_nodes = list1.size();
    
    if(interface.comportement.size()>0){
        _comportement = interface.comportement;
    }
    
    _coeffrottement = interface.coeffrottement_vec;
    _n1 = interface.side[0].neq[list1];
    _F1 = side1.F[list1];
    _F2 = side2.F[list2];
    _Wp1 = side1.Wp[list1];
    _Wp2 = side2.Wp[list2];
    _old_Wchap1 = old_side1.Wchap[list1];
    _old_Wchap2 = old_side2.Wchap[list2];
    _Fchap1.resize(nb_nodes);
    _Fchap2.resize(nb_nodes);
    _Wpchap1.resize(nb_nodes);
    _Wpchap2.resize(nb_nodes);
    
    if(interface.comp == Interface::comp_cassable_parfait or
       interface.comp == Interface::comp_cassable_elastique or
       interface.comp == Interface::comp_contact_parfait or
       interface.comp == Interface::comp_contact_elastique or
       interface.comp == Interface::comp_cohesive){
        _Wchap1.resize(nb_nodes,0.0);
        _Wchap2.resize(nb_nodes,0.0);
    }
    
    _Ep_imposee = interface.Ep_imposee[list1];
    _old_Ep_imposee = interface.old_Ep_imposee[list1];
    _Precharge = interface.precharge[list1];
    _old_Precharge = interface.old_precharge[list1];
    
    k1.kn = interface.side[0].kn;
    k1.kt = interface.side[0].kt;
    k2.kn = interface.side[1].kn;
    k2.kt = interface.side[1].kt;
    h1.kn = interface.side[0].hn;
    h1.kt = interface.side[0].ht;
    h2.kn = interface.side[1].hn;
    h2.kt = interface.side[1].ht;
    
    /// Grandeurs pour l'elasticite
    if(matprop != 0 and (interface.comp.find("Elastique") < interface.comp.size())){
        K.kn = interface.matprop->Kn;
        K.kt = interface.matprop->Kt;
        _Ep_elastique = interface.Ep_elastique;
        _old_Ep_elastique = interface.old_Ep_elastique;
    }
    
    /// Grandeurs pour l'endommagement
    if(matprop != 0 and matprop->degradable){
        _d = side1.d;
        _old_d = old_side1.d;
    }
}

void Interface::NodalState::set_node(unsigned i_node){
    id = i_node;
    LMT::Vec<unsigned> rep = LMT::range(id*DIM,(id+1)*DIM);
    
    if(interface.comportement.size()>0){
        comportement = interface.comportement[id];
    }
    
    coeffrottement = _coeffrottement[id];
    Ep_imposee = _Ep_imposee[rep];
    old_Ep_imposee = _old_Ep_imposee[rep];
    Precharge = _Precharge[rep];
    old_Precharge = _old_Precharge[rep];
    
    n1 = _n1[rep];
    F1 = _F1[rep];
    F2 = _F2[rep];
    Wp1 = _Wp1[rep];
    Wp2 = _Wp2[rep];
    old_Wchap1 = _old_Wchap1[rep];
    old_Wchap2 = _old_Wchap2[rep];
    
    Fchap1 = _Fchap1[rep];
    Fchap2 = _Fchap2[rep];
    Wpchap1 = _Wpchap1[rep];
    Wpchap2 = _Wpchap2[rep];
    
    /// Grandeurs pour l'elasticite
    if(matprop != 0 and (interface.comp.find("Elastique") < interface.comp.size())){
        K.n = n1;
        Ep_elastique = interface.Ep_elastique[rep];
        old_Ep_elastique = interface.old_Ep_elastique[rep];
    }
    
    /// Grandeurs pour l'endommagement
    if(matprop != 0 and matprop->degradable){
        d = _d[id];
        old_d = _d[id];
    }
    
    /// Maj des operateurs locaux de ddr
    k1.n = n1;
    k2.n = n1;
    h1.n = n1;
    h2.n = n1;
}

void Interface::NodalState::store_results(){
    LMT::Vec<unsigned> rep = LMT::range(id*DIM,(id+1)*DIM);
    _Fchap1[rep] = Fchap1;
    _Fchap2[rep] = Fchap2;
    _Wpchap1[rep] = Wpchap1;
    _Wpchap2[rep] = Wpchap2;
    
    if(interface.comp == Interface::comp_cassable_parfait or
       interface.comp == Interface::comp_cassable_elastique or
       interface.comp == Interface::comp_contact_parfait or
       interface.comp == Interface::comp_contact_elastique or
       interface.comp == Interface::comp_cohesive){
        _Wchap1[rep] = Wchap1;
        _Wchap2[rep] = Wchap2;
    }
    if(matprop != 0 and matprop->degradable){
        _d[id] = d;
    }
    if(matprop != 0 and (interface.comp.find("Elastique") < interface.comp.size())){
        _Ep_elastique[rep] = Ep_elastique;
    }
    if(interface.comportement.size()>0){
        _comportement[id] = comportement;
    }
}



void Interface::NodalState::save_results(){
    Side::Time& side1 = interface.side[0].t[i_time];
    Side::Time& side2 = interface.side[1].t[i_time];
    Side::Time& old_side1 = interface.side[0].t[i_time-1];
    Side::Time& old_side2 = interface.side[1].t[i_time-1];
    LMT::Vec<unsigned> &list1 = interface.side[0].ddlcorresp;
    LMT::Vec<unsigned> &list2 = interface.side[1].ddlcorresp;
    
    if(interface.comportement.size()>0){
        interface.comportement = _comportement;
    }
    
    side1.Fchap[list1] = _Fchap1;
    side2.Fchap[list2] = _Fchap2;
    side1.Wpchap[list1] = _Wpchap1;
    side2.Wpchap[list2] = _Wpchap2;
    
    if(interface.comp == Interface::comp_cassable_parfait or
       interface.comp == Interface::comp_cassable_elastique or
       interface.comp == Interface::comp_contact_parfait or
       interface.comp == Interface::comp_contact_elastique or
       interface.comp == Interface::comp_cohesive){
        side1.Wchap[list1] = _Wchap1;
        side2.Wchap[list2] = _Wchap2;
    }
    
    /// Grandeurs pour l'elasticite
    if(matprop != 0 and (interface.comp.find("Elastique") < interface.comp.size())){
        interface.Ep_elastique = _Ep_elastique;
    }
    
    /// Grandeurs pour l'endommagement
    if(matprop != 0 and matprop->degradable){
        side1.d = _d;
    }
}
