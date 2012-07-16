#ifndef THERMALLOAD_H
#define THERMALLOAD_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"
#include "../COMPUTE/DataUser.h"
#include "Sst.h"
#include "../../LMT/include/mesh/tetra.h"

/** Cette classe est le pendant de la classe Boundary (c.f. DEFINITION/Boundary.h)
 * Une seule instance de cette classe suffit :
 * elle reprend l'expression de la temperature imposée
 */
struct ThermalLoad{
    struct AssigneOnElement {
        template<class TE>
        void operator() (TE &e, ThermalLoad &Tload,Codegen::Ex::MapExNum &values) const {
            Point G = LMT::center(e);
            values[Tload.x.self_ex] = G[0];
            values[Tload.y.self_ex] = G[1];
            #if DIM == 3
            values[Tload.z.self_ex] = G[2];
            #endif
            e.deltaT = Tload.deltaT.updateValue(values);// * LMT::measure(e);// Inutile on dirait...
        }
    };
    
    ThermalLoad();
    
    void read_data_user(const DataUser &data_user);
    void prepareParameters();
    void updateParameters();
    
    void apply_on_sst(Sst &S){
        Codegen::Ex::MapExNum values = parameters.getParentsValues();
        LMT::apply(S.mesh->elem_list,AssigneOnElement(),*this, values);
    }
    
    void affiche() const;
    
    ParameterGroup parameters;
    MainParameter x;
    MainParameter y;
    MainParameter z;
    UserParameter deltaT;
};

#endif //THERMALLOAD_H