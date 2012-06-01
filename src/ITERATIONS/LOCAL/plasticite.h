#ifndef PLASTICITE_H
#define PLASTICITE_H

#include "../../MAILLAGE/elements_variables_accessor.h"
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"
#include "../../DEFINITIONS/Interface.h"
#include "../../UTILS/utils_2.h"
//#include <boost/bind/bind_template.hpp>

void calcul_plasticite(Sst &S, Process &process);

#endif //PLASTICITE_H