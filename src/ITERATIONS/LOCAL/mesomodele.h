#ifndef MESOMODELE_H
#define MESOMODELE_H

#include "../../DEFINITIONS/Process.h"
#include "../../DEFINITIONS/structure_typedef.h"


struct reactualisation_rigidite {
    void operator()(Sst &S, VecInterfaces &Inter, const Process &process, const DataUser &data_user) const;
};

void calcul_mesomodele(Sst &S, Process &process);

#endif //MESOMODELE_H
