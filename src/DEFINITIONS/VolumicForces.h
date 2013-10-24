#ifndef VOLUMICFORCES_H
#define VOLUMICFORCES_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"
#include "../COMPUTE/DataUser.h"
#include "Sst.h"
#include "../../LMT/include/mesh/tetra.h"

/** Cette classe est le pendant de la classe Boundary (c.f. DEFINITION/Boundary.h)
 * Une seule instance de cette classe suffit :
 * toutes les expressions d'efforts volumiques y sont concatenees.
 */
struct VolumicForces{
    struct AssignOnElement {
        template<class TE>
        void operator() (TE &e, Scalar density, VolumicForces &Fvol,Codegen::Ex::MapExNum &values) const {
            Point G = LMT::center(e);
            values[Fvol.x.self_ex] = G[0];
            values[Fvol.y.self_ex] = G[1];
            #if DIM == 3
            values[Fvol.z.self_ex] = G[2];
            #endif
            for(unsigned d=0;d<DIM;++d){
                //e.f_vol_e[d] = density * Fvol.force[d].updateValue(values);// * LMT::measure(e);// Inutile on dirait...
                e.f_vol_e[d] = Fvol.force[d].updateValue(values);// * LMT::measure(e);// Inutile on dirait...
                //std::cout << density << endl;
                //std::cout << e.f_vol_e[d] << endl;
            }
        }
    };
    
    VolumicForces();
    
    void read_data_user(const DataUser &data_user);
    void prepareParameters();
    void updateParameters();
    
    void apply_on_sst(Sst &S){
        
        std::cout << "Application des efforts volumiques sur S (id " << S.id << "):" << std::endl;
        for(int i = 0; i < DIM; i++){
            std::cout << "\t" << force[i].str_expr << std::endl;
        }
        //
        Codegen::Ex::MapExNum values = parameters.getParentsValues();
        LMT::apply(S.mesh->elem_list,AssignOnElement(),S.mesh->density,*this, values);
    }
    
    void affiche() const;
    
    ParameterGroup parameters;
    MainParameter x;
    MainParameter y;
    MainParameter z;
    LMT::Vec<UserParameter,3> force;
};

#endif //VOLUMICFORCES_H