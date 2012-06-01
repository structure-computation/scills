#include "elements_variables_accessor.h"


void upload_q(Sst &S,Vec<TYPEREEL> &sto){
    S.f->get_result() = sto;
}

void upload_q(Sst &S,Sst::Time &t){
    S.f->get_result() = t.q;
}


void upload_q(Sst &S,unsigned pt){
    S.f->get_result() = S.t[pt].q;
}



//*
///---------------------------- SIGMA ----------------------------------
/// Recuperation depuis le maillage
void download_sigma(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_sigma(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_sigma(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_sigma(),sto,i_elem);
}
//*/
/*
///---------------------------- SIGMA VON MISES ----------------------------
/// Recuperation depuis le maillage
void download_sigma_von_mises(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_sigma_von_mises(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_sigma_von_mises(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_sigma_von_mises(),sto,i_elem);
}
//*/
//*
///---------------------------- EPSILON ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_epsilon(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_epsilon(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_epsilon(),sto,i_elem);
}
//*/
/*
///---------------------------- EPSILON E ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon_e(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_epsilon_e(),sto,i_elem);
}

/// Envoi vers le maillage
void upload_epsilon_e(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_epsilon_e(),sto,i_elem);
}
//*/

///---------------------------- EPSILON P ---------------------------------
/// Recuperation depuis le maillage
void download_epsilon_p(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_epsilon_p(),sto,i_elem);
}

void download_epsilon_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_epsilon_p(),t.epsilon_p,i_elem);
}

void download_epsilon_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_epsilon_p(),S.t[pt].epsilon_p,i_elem);
}

/// Envoi vers le maillage
void upload_epsilon_p(Sst &S,Vec<Vec<TYPEREEL,DIM*(DIM+1)/2> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_epsilon_p(),sto,i_elem);
}

void upload_epsilon_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_epsilon_p(),t.epsilon_p,i_elem);
}

void upload_epsilon_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_epsilon_p(),S.t[pt].epsilon_p,i_elem);
}

//*
///---------------------------     P     ----------------------------------
/// Recuperation depuis le maillage
void download_p(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_p(),sto,i_elem);
}

void download_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_p(),t.p,i_elem);
}

void download_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_p(),S.t[pt].p,i_elem);
}

/// Envoi vers le maillage
void upload_p(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_p(),sto,i_elem);
}

void upload_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_p(),t.p,i_elem);
}

void upload_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_p(),S.t[pt].p,i_elem);
}
//*/
//*
///-------------------------------  R_p  -------------------------------------------
/// Recuperation depuis le maillage
void download_R_p(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_R_p(),sto,i_elem);
}

void download_R_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_R_p(),t.R_p,i_elem);
}

void download_R_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_R_p(),S.t[pt].R_p,i_elem);
}

/// Envoi vers le maillage
void upload_R_p(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_R_p(),sto,i_elem);
}

void upload_R_p(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_R_p(),t.R_p,i_elem);
}

void upload_R_p(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_R_p(),S.t[pt].R_p,i_elem);
}
//*/
//*
///-------------------------------   d1   -------------------------------------------
/// Recuperation depuis le maillage
void download_d1(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d1(),sto,i_elem);
}

void download_d1(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d1(),t.d1,i_elem);
}

void download_d1(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d1(),S.t[pt].d1,i_elem);
}

/// Envoi vers le maillage
void upload_d1(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d1(),sto,i_elem);
}

void upload_d1(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d1(),t.d1,i_elem);
}

void upload_d1(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d1(),S.t[pt].d1,i_elem);
}
//*/
//*
///-------------------------------   d2   -------------------------------------------
/// Recuperation depuis le maillage
void download_d2(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d2(),sto,i_elem);
}

void download_d2(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d2(),t.d1,i_elem);
}

void download_d2(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_d2(),S.t[pt].d1,i_elem);
}

/// Envoi vers le maillage
void upload_d2(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d2(),sto,i_elem);
}

void upload_d2(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d2(),t.d1,i_elem);
}

void upload_d2(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_d2(),S.t[pt].d1,i_elem);
}
//*/
//*
///-------------------------------   df   -------------------------------------------
/// Recuperation depuis le maillage
void download_df(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_df(),sto,i_elem);
}

void download_df(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_df(),t.d1,i_elem);
}

void download_df(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_df(),S.t[pt].d1,i_elem);
}

/// Envoi vers le maillage
void upload_df(Sst &S,Vec<TYPEREEL> &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_df(),sto,i_elem);
}

void upload_df(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_df(),t.d1,i_elem);
}

void upload_df(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_df(),S.t[pt].d1,i_elem);
}
//*/
//*
///-------------------------------   Yd   -------------------------------------------
/// Recuperation depuis le maillage
void download_Yd(Sst &S,Vec<Vec<TYPEREEL> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_Yd(),sto,i_elem);
}

void download_Yd(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_Yd(),t.Yd,i_elem);
}

void download_Yd(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Download_Yd(),S.t[pt].Yd,i_elem);
}

/// Envoi vers le maillage
void upload_Yd(Sst &S,Vec<Vec<TYPEREEL> > &sto){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_Yd(),sto,i_elem);
}

void upload_Yd(Sst &S,Sst::Time &t){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_Yd(),t.Yd,i_elem);
}

void upload_Yd(Sst &S,unsigned pt){
    unsigned i_elem = 0;
    apply(S.mesh->elem_list,__Upload_Yd(),S.t[pt].Yd,i_elem);
}
//*/
