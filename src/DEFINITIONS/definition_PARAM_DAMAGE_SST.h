#ifndef PARAM_DAMAGE_SST_H
#define PARAM_DAMAGE_SST_H
#include "mesh/mesh.h"

#include "meshcaracinter.h" //utilise pour la class PARAM_DAMAGE_SST (uniquement pour le mesomodele et les sst)

template<unsigned dim_, class T_>
struct PARAM_DAMAGE_SST {
    void free() {
        mesh.free();
        d.free();
        dmax.free();
        Y.free();
        Ymax.free();
    }
    typedef T_ T;
    static const unsigned dim=dim_;
    typedef Mesh<Meshcaracinter<dim,dim-1> > TMESH;
    PARAM_DAMAGE_SST() {
        Yo = 0.0961;
        Yop = 0.0961;
        Ysp = 0.70;
        Yc = 10.8;
        Ycp = 10.8;
        b = 2.5;
        d.set(0);
        dmax.set(0);
        Y.set(0);
        Ymax.set(0);
    }
    T Yo,Yop, Ysp, Yc, Ycp, b, epaisseur;
    Vec<T,2> d, dmax, Y, Ymax;
    TMESH mesh;
};


#endif //PARAM_DAMAGE_SST_H


