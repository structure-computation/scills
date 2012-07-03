#ifndef VOLUMICFORCES_H
#define VOLUMICFORCES_H

#include "main_typedef.h"
#include "../PARAMETERS/Parameters.h"
#include "../COMPUTE/DataUser.h"


/** Cette classe est le pendant de la classe Boundary (c.f. DEFINITION/Boundary.h)
 * Une seule instance de cette classe suffit :
 * toutes les expressions d'efforts volumiques y sont concatenees.
 */
struct VolumicForces{
    VolumicForces();
    
    void read_data_user(const DataUser &data_user);
    void prepareParameters();
    void updateParameters();
    void affiche() const;
    
    ParameterGroup parameters;
    MainParameter x;
    MainParameter y;
    MainParameter z;
    LMT::Vec<UserParameter,3> force;
};

#endif //VOLUMICFORCES_H