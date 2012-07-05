#include "Process.h"
#include "Sst.h"
#include "SstCarac_InterCarac.h"
#include "../UTILS/utils_2.h"

/// constructeur de la formulation pour la sous-structure
Sst::Sst() : 
pb(*mesh.m,true),
plastique(false),
endommageable(false),
matprop(0),
update_operator(true)
{}

Sst::~Sst() {free();}    ///destructeur de la SST

void Sst::free(){
    vois.free();
    edge.free();
    LE.free();
    /// K est liberee apres factorisation de toute facon et non calculer où y a pas besoin
    fvol.free();
    t.free();
    /// matprop est libere par le Process
}

/*void Sst::resize_results(unsigned nb_t,bool allocate_t,unsigned nb_t_post,bool allocate_t_post,bool recopie_t_post){
    t.resize(nb_t);
    if(recopie_t_post){
        t = t_post;
    }
    if(allocate_t){
        for(unsigned i = 0; i < nb_t; i++){
            t[i].allocate(*this);
        }
    }
}*/

void Sst::read_data_user(int index,const DataUser &data_user,GeometryUser &geometry_user){/*
    const DataUser::Json_pieces &sst = data_user.pieces_vec[index];
    id = sst.id_in_calcul;
    num = sst.id_in_calcul;
    id_material = data_user.find_materials_index(sst.material_id);
    
    /// Chargement du maillage associe
    mesh.name=sst.name;
    mesh.load(geometry_user, id);
    mesh.load();
    mesh->update_skin();
    geometry_user.find_group_elements(id)->nb_elements_skin=mesh->skin.elem_list.size();   /// A REVOIR
    geometry_user.find_group_elements(id)->nb_nodes_skin=mesh->skin.node_list.size();      /// A REVOIR
    
    int nb_nodes_by_element_skin=(geometry_user.patterns.find_type(geometry_user.find_group_elements(id)->pattern_base_id)).nb_nodes_by_sides;
    geometry_user.find_group_elements(id)->local_connectivities_skin.resize(nb_nodes_by_element_skin);
    for(int ne=0;ne<nb_nodes_by_element_skin;ne++){
        geometry_user.find_group_elements(id)->local_connectivities_skin[ne].resize(mesh->skin.elem_list.size());
    }
    apply(mesh->skin.elem_list,ConvertMeshConnectivitiesSkin(),geometry_user, id);
    mesh.unload();//*/
}

/// Trouver une sst à partir de son id----------------------------------------------
Sst* Sst::find_sst(LMT::Vec<Sst> &S,int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return &S[i_group];
            break;
        }
    }
    std::cerr << "Impossible de trouver la Sst " << id_ << std::endl;
    assert(0);
}

/// Trouver l'index d'une sst à partir de son id -----------------------------------
int Sst::find_index_sst(LMT::Vec<Sst> &S, int id_) {
    for (int i_group=0; i_group<S.size(); i_group++) {
        if (S[i_group].id == id_) {
            return i_group;
            break;
        }
    }
    std::cerr << "Impossible de trouver la Sst " << id_ << std::endl;
    assert(0);
}

void Sst::apply_behavior(){
    mesh->density    = matprop->density    ;
    mesh->deltaT     = matprop->deltaT     ;
    mesh->resolution = matprop->resolution ;
    
    if (matprop->type.find("isotrope")<matprop->type.size()) {
        ///formulation isotrope
        mesh->elastic_modulus = matprop->elastic_modulus ;
        mesh->poisson_ratio   = matprop->poisson_ratio   ;
        mesh->alpha           = matprop->alpha           ;
    } else if (matprop->type.find("orthotrope")<matprop->type.size()) {
        ///formulation orthotrope
        mesh->elastic_modulus_1 = matprop->elastic_modulus_1;
        mesh->elastic_modulus_2 = matprop->elastic_modulus_2;
        mesh->elastic_modulus_3 = matprop->elastic_modulus_3;
        mesh->poisson_ratio_12  = matprop->poisson_ratio_12 ;
        mesh->poisson_ratio_13  = matprop->poisson_ratio_13 ;
        mesh->poisson_ratio_23  = matprop->poisson_ratio_23 ;
        mesh->shear_modulus_12  = matprop->shear_modulus_12 ;
        mesh->shear_modulus_13  = matprop->shear_modulus_13 ;
        mesh->shear_modulus_23  = matprop->shear_modulus_23 ;
        mesh->v1                = matprop->v1               ;
        mesh->v2                = matprop->v2               ;
        mesh->deltaT            = matprop->deltaT           ;
        mesh->resolution        = matprop->resolution       ;
        mesh->alpha_1           = matprop->alpha_1          ;
        mesh->alpha_2           = matprop->alpha_2          ;
        mesh->alpha_3           = matprop->alpha_3          ;
        mesh->f_vol             = matprop->f_vol            ;
    }
    
    if(matprop->comp.find("pl")<matprop->comp.size()){
        mesh->plast_ecrouissage_init = matprop->plast_ecrouissage_init;
        mesh->plast_ecrouissage_mult = matprop->plast_ecrouissage_mult;
        mesh->plast_ecrouissage_expo = matprop->plast_ecrouissage_expo;
        mesh->plast_cinematique_coef = matprop->plast_cinematique_coef;
    }
    
    if(matprop->comp.find("en")<matprop->comp.size()){
        mesh->Yo           = matprop->Yo;
        mesh->b_c          = matprop->b_c;
    }
    
    if(matprop->comp.find("mesomodele")<matprop->comp.size()){
        mesh->Yo           = matprop->Yo;
        mesh->Yc           = matprop->Yc;
        mesh->Ycf          = matprop->Ycf;
        mesh->dmax         = matprop->dmax;
        mesh->b_c          = matprop->b_c;
        mesh->effet_retard = matprop->effet_retard;
        mesh->a            = matprop->a;
        mesh->tau_c        = matprop->tau_c;
    }
}



Sst::Time::Time(){
}

Sst::Time::~Time(){
    free();
}

void Sst::Time::allocate(Sst &S){
    unsigned nbddl = S.mesh.node_list_size*DIM;
    unsigned nbelem = S.mesh.elem_list_size;
    /// Deplacements
    q.resize(nbddl);
    q.set(0.0);
    /// Comportement plastique
    if(S.f == S.pb.formulation_plasticity_isotropy_stat_Qstat){
        /// Plasticite cumulee
        p.resize(nbelem);
        p.set(0.0);
        /// Ecrouissage isotrope
        R_p.resize(nbelem);
        R_p.set(0.0);
        /// Deformations plastiques
        epsilon_p.resize(nbelem);
        for( unsigned i = 0; i < nbelem; ++i ) {
            epsilon_p[i].set(0.0);
        }
        /// Origines des domaines elastiques
        if(S.matprop->type_plast.find("cinematique") < S.matprop->type_plast.size()){
            X_p.resize(nbelem);
            for( unsigned i = 0; i < nbelem; ++i ) {
                X_p[i].set(0.0);
            }
        }
    }
    /// Comportement endommageable
    if(S.f == S.pb.formulation_elasticity_damageable_isotropy_stat_Qstat or S.f == S.pb.formulation_mesomodele){
        /// Endommagement standart ou micro-fissuration du mesomodele
        d1.resize(nbelem);
        d1.set(0.0);
    }
    /// Comportement mesomodele
    if(S.f == S.pb.formulation_mesomodele){
        /// Decohesion fibres/matrice du mesomodele
        d2.resize(nbelem);
        d2.set(0.0);
        /// Rupture fragile des fibres du mesomodele
        df.resize(nbelem);
        df.set(0.0);
        /// Forces associees
        Yd.resize(nbelem);
        for( unsigned i = 0; i < nbelem; ++i )
            Yd[i].set(0.0);
    }
}

void Sst::Time::free(){
    if(q.size()!=0){q.free();}
    if(p.size()!=0){p.free();}
    if(R_p.size()!=0){R_p.free();}
    if(epsilon_p.size()!=0){epsilon_p.free();}
    if(X_p.size()!=0){X_p.free();}
    if(d1.size()!=0){d1.free();}
    if(d2.size()!=0){d2.free();}
    if(df.size()!=0){df.free();}
    if(Yd.size()!=0){Yd.free();}
}

void Sst::Time::affiche() const {
    std::cout << "********************** Debug Sst::Time : ********************" << std::endl;
    std::cout << "q         : " << q.size()         << std::endl;
    std::cout << "p         : " << p.size()         << std::endl;
    std::cout << "R_p       : " << R_p.size()       << std::endl;
    std::cout << "epsilon_p : " << epsilon_p.size() << std::endl;
    std::cout << "d1        : " << d1.size()        << std::endl;
    std::cout << "d2        : " << d2.size()        << std::endl;
    std::cout << "df        : " << df.size()        << std::endl;
    std::cout << "Yd        : " << Yd.size()        << std::endl;
    std::cout << "*************************************************************" << std::endl;
    std::cout << std::endl;
}

void Sst::calc_SST_Correspddl() {
    for(unsigned j=0;j<edge.size();++j) {
        edge[j].repddledge.resize(edge[j].mesh->node_list.size()*DIM);
        IntersectionCarac<EdgeMesh::TNodeList, SstMesh::TM::TNodeList >::T inter=
        intersection_ptr(edge[j].mesh->node_list,mesh->node_list,Function<DistBetweenNodes>() < 0.0001);
        
        Vec<unsigned> repnode;
        repnode.resize(edge[j].mesh->node_list.size());
        
        if (inter.size() != repnode.size()) {
            std::cerr << "Calc_SST_Correspddl - attention noeud non repere - probleme de correspondance - cote " << j << endl;
            std::cerr << "meshini - size node " << mesh.node_list_size << " numero " << num << endl;
            std::cerr << "intersize "<< inter.size() << " edge size "  << edge[j].mesh->node_list.size()<< " internum " << edge[j].internum << endl;
            assert(0);
        }
        
        for(unsigned i=0;i<inter.size();++i)
            repnode[inter[i].first->number_in_original_mesh()]=inter[i].second->number_in_original_mesh();
        for(unsigned i=0;i<repnode.size();++i)
            edge[j].repddledge[range(i*DIM,(i+1)*DIM)]=range(repnode[i]*DIM,(repnode[i]+1)*DIM);
    }
}

void Sst::affiche() const {
    std::cout << "************************* Debug Sst : ***********************" << std::endl;
    debug("G              ",G);
    debug("measure        ",measure);
    debug("num            ",num);
    debug("id             ",id);
    debug("num_proc       ",num_proc);
    debug("vois           ",vois);
    debug("id_material    ",id_material);
    debug("typmat         ",typmat);
    debug("plastique      ",plastique);
    debug("endommageable  ",endommageable);
    debug("update_operator",update_operator);
    debug("nb_macro       ",nb_macro);
    debug("nb_macro_espace",nb_macro_espace);
    debug("nb nodes       ",mesh.node_list_size);
    debug("nb elements    ",mesh.elem_list_size);
    //std::cout << "formulation : " << f << std::endl;
    std::cout << "edge : " << edge.size() << std::endl;
    for(int i = 0; i < edge.size(); i++){
        edge[i].affiche();
    }
    std::cout << "t : " << t.size() << std::endl;
    for(int i = 0; i < t.size(); i++){
        t[i].affiche();
    }
    std::cout << "t_post : " << t.size() << std::endl;
    for(int i = 0; i < t.size(); i++){
        t[i].affiche();
    }
    std::cout << "*************************************************************" << std::endl;
}
