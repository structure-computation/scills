#ifndef MESOMODELE_H
#define MESOMODELE_H

#include "../../DEFINITIONS/Interface.h"
#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/Sst.h"


struct reactualisation_rigidite {
    void operator()(Sst &S, Vec<Interface> &Inter, const Process &process, const DataUser &data_user) const;
};

void calcul_mesomodele(Sst &S, Process &process);

#endif //MESOMODELE_H
