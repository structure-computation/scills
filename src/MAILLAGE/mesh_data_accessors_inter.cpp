#include "mesh_data_accessors_inter.h"


///------------------------------------ F ----------------------------------------
/// Recuperation depuis le maillage
void download_F(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_F(),sto,i_elem);
}

void download_F(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_F(),t.F,i_elem);
}

void download_F(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_F(),side.t[pt].F,i_elem);
}

/// Envoi vers le maillage
void upload_F(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_F(),sto,i_elem);
}

void upload_F(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_F(),t.F,i_elem);
}

void upload_F(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_F(),side.t[pt].F,i_elem);
}


///------------------------------------ W ----------------------------------------
/// Recuperation depuis le maillage
void download_W(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_W(),sto,i_elem);
}

void download_W(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_W(),t.W,i_elem);
}

void download_W(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_W(),side.t[pt].W,i_elem);
}

/// Envoi vers le maillage
void upload_W(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_W(),sto,i_elem);
}

void upload_W(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_W(),t.W,i_elem);
}

void upload_W(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_W(),side.t[pt].W,i_elem);
}


///----------------------------------- dWn ---------------------------------------
/// Recuperation depuis le maillage
void download_dWn(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_dWn(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_dWn(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_dWn(),sto,i_elem);
}

///----------------------------------- dWt ---------------------------------------
/// Recuperation depuis le maillage
void download_dWt(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_dWt(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_dWt(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_dWt(),sto,i_elem);
}


///------------------------------- dissipation -----------------------------------
/// Recuperation depuis le maillage
void download_dissipation(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_dissipation(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_dissipation(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_dissipation(),sto,i_elem);
}


///------------------------------------ d ----------------------------------------
/// Recuperation depuis le maillage
void download_d(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_d(),sto,i_elem);
}

void download_d(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_d(),t.d,i_elem);
}

void download_d(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Download_d(),side.t[pt].d,i_elem);
}

/// Envoi vers le maillage
void upload_d(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_d(),sto,i_elem);
}

void upload_d(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_d(),t.d,i_elem);
}

void upload_d(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(side.mesh->elem_list,__Upload_d(),side.t[pt].d,i_elem);
}

