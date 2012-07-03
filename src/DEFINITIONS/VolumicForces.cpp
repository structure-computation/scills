#include "VolumicForces.h"

VolumicForces::VolumicForces():
parameters(),
x("x"),
force(Vec<UserParameter,3>(UserParameter("gamma_x"),UserParameter("gamma_y"),UserParameter("gamma_z"))),
y("y"),
z("z")
{
    parameters.addParameter(&x);
    parameters.addParameter(&y);
    parameters.addParameter(&z);
    for(int d = 0; d < DIM; d++){
        parameters.addParameter(&(force[d]));
    }
}

void VolumicForces::read_data_user(const DataUser &data_user){
    /// On concatene les expressions sur chaque direction:((gamma1)*(dx1))+((gamma1)*(dx2))++...
    /// En cas de - uniaire, on rajoute un 0 apres la parenthese ouvrante : ...+((0-10)*(...))+...
    Sc2String F[DIM];
    for(int d = 0; d < DIM; d++){
        F[d] << "((0";
    }
    /// Au cas ou aucun effort volumique ne serait defini
    if(data_user.volumic_forces_vec.size() == 0){
        for(int d = 0; d < DIM; d++){
            force[d].setExpression("0");
        }
        return;
    }
    for(unsigned i = 0; i < data_user.volumic_forces_vec.size(); i++){
        /// Preparation pour un nouvel effort
        for(int d = 0; d < DIM; d++){
            F[d] << "))+((";
        }
        /// Insertion de gamma
        const Sc2String &g = data_user.volumic_forces_vec[i].gamma;
        if(g[0] == '-'){
            for(int d = 0; d < DIM; d++){
                F[d] << "0";
            }
        }
        for(int d = 0; d < DIM; d++){
            F[d] << g << ")*(";
        }
        /// Insertion de dx, dy, dz
        const Sc2String &fx = data_user.volumic_forces_vec[i].dx;
        if(fx[0] == '-'){
            F[0] << "0";
        }
        F[0] << fx;
        const Sc2String &fy = data_user.volumic_forces_vec[i].dy;
        if(fy[0] == '-'){
            F[1] << "0";
        }
        F[1] << fy;
        #if DIM == 3
        const Sc2String &fz = data_user.volumic_forces_vec[i].dz;
        if(fz[0] == '-'){
            F[2] << "0";
        }
        F[2] << fz;
        #endif
    }
    for(int d = 0; d < DIM; d++){
        F[d] << "))";                   /// Finalisation de l'expression
        std::cout << F[d] << std::endl;
        force[d].setExpression(F[d]);  /// Assignation dans force
    }
}

void VolumicForces::prepareParameters(){
    parameters.prepareParameters();
}

void VolumicForces::updateParameters(){
    /// Remettre ici le code de maj des CL
}

void VolumicForces::affiche() const{
    std::cout << "*************************************************************" << std::endl;
    std::cout << "******************* Debug VolumicForces : *******************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
    LMT::Vec<Sc2String,3> dir;
    dir[0] = "x";
    dir[1] = "y";
    dir[2] = "z";
    for(unsigned d = 0; d < DIM; d++){
        std::cout << "    f" << dir[d] << " : " << force[d].str_expr << " = " << force[d] << std::endl;
    }
    std::cout << "*************************************************************" << std::endl;
    std::cout << "*************************************************************" << std::endl;
}