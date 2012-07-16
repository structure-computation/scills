#ifndef PLASTICITE_H
#define PLASTICITE_H

#include "../../MAILLAGE/mesh_data_accessors_sst.h"
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILS/utils_2.h"
//#include <boost/bind/bind_template.hpp>

void calcul_plasticite(Sst &S, Process &process);

#endif //PLASTICITE_H
