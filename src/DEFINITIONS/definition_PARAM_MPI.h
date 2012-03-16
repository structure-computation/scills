#ifndef PARAM_MULTI_MPI_H
#define PARAM_MULTI_MPI_H
using namespace LMT;

/** \ingroup Parametres
\brief Parmamètres associés aux données MPI
*/
struct MULTI_MPI{
  ///Constructeur : valeurs par défaut pour les paramètres
  MULTI_MPI(){
//    rank=0;
//    size=0;
  }
//int rank;
//int size;
Vec<Vec<int> > repartition_sst;
Vec<unsigned> listsst;
Vec<unsigned> listinter;


struct INTERTOEXCHANGE1 {
int num;//numero de l interface a envoyer
int sidetosend;//numero du cote a envoyer
int to;//numero du processeur a qui on envoie
};
Vec<INTERTOEXCHANGE1> intertoexchange;
struct INTERTOEXCHANGEBYPRO1 {
int to;
Vec<INTERTOEXCHANGE1> inter;
};
Vec<INTERTOEXCHANGEBYPRO1> intertoexchangebypro;

Vec<INTERTOEXCHANGE1> intertoexchangeformaster;

typedef INTERTOEXCHANGE1 INTERTOEXCHANGE;
typedef INTERTOEXCHANGEBYPRO1 INTERTOEXCHANGEBYPRO;

};

#endif //PARAM_MULTI_MPI_H

