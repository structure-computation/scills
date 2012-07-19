#ifndef MESH_DATA_ACCESSORS_INTER_H
#define MESH_DATA_ACCESSORS_INTER_H

#include "../DEFINITIONS/main_typedef.h"
#include "../DEFINITIONS/Interface.h"


///------------------------------------ F ----------------------------------------
/// Recuperation depuis le maillage
void download_F(Interface::Side &side,Vector &sto);
void download_F(Interface::Side &side,Interface::Side::Time &t);
void download_F(Interface::Side &side,unsigned pt);

struct __Download_F{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.F[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_F(Interface::Side &side,Vector &sto);
void upload_F(Interface::Side &side,Interface::Side::Time &t);
void upload_F(Interface::Side &side,unsigned pt);

struct __Upload_F{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.F[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};


///------------------------------------ W ----------------------------------------
/// Recuperation depuis le maillage
void download_W(Interface::Side &side,Vector &sto);
void download_W(Interface::Side &side,Interface::Side::Time &t);
void download_W(Interface::Side &side,unsigned pt);

struct __Download_W{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.W[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_W(Interface::Side &side,Vector &sto);
void upload_W(Interface::Side &side,Interface::Side::Time &t);
void upload_W(Interface::Side &side,unsigned pt);

struct __Upload_W{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.W[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};


///----------------------------------- dWn ---------------------------------------
/// Recuperation depuis le maillage
void download_dWn(Interface::Side &side,Vector &sto);

struct __Download_dWn{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.dWn[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_dWn(Interface::Side &side,Vector &sto);

struct __Upload_dWn{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.dWn[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};


///----------------------------------- dWt ---------------------------------------
/// Recuperation depuis le maillage
void download_dWt(Interface::Side &side,Vector &sto);

struct __Download_dWt{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.dWt[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_dWt(Interface::Side &side,Vector &sto);

struct __Upload_dWt{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.dWt[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};


///------------------------------- dissipation -----------------------------------
/// Recuperation depuis le maillage
void download_dissipation(Interface::Side &side,Vector &sto);

struct __Download_dissipation{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.dissipation[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_dissipation(Interface::Side &side,Vector &sto);

struct __Upload_dissipation{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.dissipation[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};


///------------------------------------ d ----------------------------------------
/// Recuperation depuis le maillage
void download_d(Interface::Side &side,Vector &sto);
void download_d(Interface::Side &side,Interface::Side::Time &t);
void download_d(Interface::Side &side,unsigned pt);

struct __Download_d{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        sto[i_elem] = e.d;
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_d(Interface::Side &side,Vector &sto);
void upload_d(Interface::Side &side,Interface::Side::Time &t);
void upload_d(Interface::Side &side,unsigned pt);

struct __Upload_d{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        e.d = sto[i_elem];
        i_elem++;
    }
};


/*


///-------------------------------------------------------------------------------
/// Recuperation depuis le maillage
void download_*xxx*(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Download_*xxx*(),sto,i_elem);
}

void download_*xxx*(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Download_*xxx*(),t.*xxx*,i_elem);
}

void download_*xxx*(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Download_*xxx*(),I.t[pt].*xxx*,i_elem);
}

struct __Download_*xxx*{
    template<typename Elem,typename Stockage>
    void operator()(const Elem &e, Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            sto[i_elem*DIM+i] = e.*xxx*[i];
        }
        i_elem++;
    }
};

/// Envoi vers le maillage
void upload_*xxx*(Interface::Side &side,Vector &sto){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Upload_*xxx*(),sto,i_elem);
}

void upload_*xxx*(Interface::Side &side,Interface::Side::Time &t){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Upload_*xxx*(),t.*xxx*,i_elem);
}

void upload_*xxx*(Interface::Side &side,unsigned pt){
    unsigned i_elem = 0;
    apply(I.mesh->elem_list,__Upload_*xxx*(),I.t[pt].*xxx*,i_elem);
}

struct __Upload_*xxx*{
    template<typename Elem,typename Stockage>
    void operator()(Elem &e, const Stockage &sto, unsigned &i_elem) const {
        for(int i = 0; i < DIM; i++){
            e.*xxx*[i] = sto[i_elem*DIM+i];
        }
        i_elem++;
    }
};
//*/

#endif //MESH_DATA_ACCESSORS_INTER_H