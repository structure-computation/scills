#ifndef PARALLELISATIONDATA_H
#define PARALLELISATIONDATA_H

#include "main_typedef.h"
#include "../COMPUTE/DataUser.h"
#include "../MPI/mpi_lmt_functions.h"

/** \ingroup Parametres
\brief Parmamètres associes aux donnees MPI
*/
struct ParallelisationData{
    //attributs==============================================================================================
    
    int size;                           /// Nombre de processeurs alloues au calcul
    int rank;                           /// Numero du processeur gerant ce programme
    int nb_threads;                     /// Nombre de threads par machine (1 par defaut)
    
    Vec<Vec<int> > repartition_sst;     /// Listes des sst gerees par chaque programme MPI
    Vec<unsigned> listsst;      /// Liste des SST gerees par cette instance MPI
    Vec<unsigned> listinter;    /// Liste des interfaces gerees par cette instance MPI


    struct INTERTOEXCHANGE {
        int num;            /// Numero de l interface a envoyer
        int sidetosend;     /// Numero du cote a envoyer
        int to;             /// Numero du processeur a qui on envoie
    };
    
    struct INTERTOEXCHANGEBYPRO {
        int to;
        Vec<INTERTOEXCHANGE> inter;
    };
    
    Vec<INTERTOEXCHANGE> intertoexchange;
    Vec<INTERTOEXCHANGE> intertoexchangeformaster;
    Vec<INTERTOEXCHANGEBYPRO> intertoexchangebypro;

    //methodes===============================================================================================
    /// Constructeur (par defaut)
    ParallelisationData():size(1),rank(0),nb_threads(1){}
    void read_data_user(const Metil::DataUser &data_user);
    /// Synchronise les differentes instances MPI du programme (si necessaire)
    void synchronisation(){if(size>1) {MPI_Barrier(MPI_COMM_WORLD);}}
    /// Indique si cette instance du programme gere des problemes locaux
    bool is_local_cpu() const {return (rank>0 or size==1);}
    /// Indique si cette instance du programme est le master
    bool is_master_cpu() const {return rank==0;}
    /// Indique si plusieurs instances du programme sont en cours d'execution
    bool is_multi_cpu() const {return size>1;}
};

#endif //PARALLELISATIONDATA_H

