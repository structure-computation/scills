#ifndef PARAM_MULTI_MPI_H
#define PARAM_MULTI_MPI_H

#include "../../LMT/include/containers/vec.h"
using namespace LMT;

/** \ingroup Parametres
\brief Parmamètres associés aux données MPI
*/
struct MULTI_MPI{
    Vec<Vec<int> > repartition_sst;
    Vec<unsigned> listsst;
    Vec<unsigned> listinter;


    struct INTERTOEXCHANGE {
        int num;            ///numero de l interface a envoyer
        int sidetosend;     ///numero du cote a envoyer
        int to;             ///numero du processeur a qui on envoie
    };
    
    struct INTERTOEXCHANGEBYPRO {
        int to;
        Vec<INTERTOEXCHANGE> inter;
    };
    
    Vec<INTERTOEXCHANGE> intertoexchange;
    Vec<INTERTOEXCHANGEBYPRO> intertoexchangebypro;
    Vec<INTERTOEXCHANGE> intertoexchangeformaster;

};

#endif //PARAM_MULTI_MPI_H

