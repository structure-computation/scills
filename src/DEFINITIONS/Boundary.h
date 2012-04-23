#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "../COMPUTE/DataUser.h"
#include "../UTILS/Sc2String.h"
#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/codegen/codegen.h"
#include "../../LMT/include/containers/basicops.h"
#include "Process.h"
//#include "TimeParameters.h"

/** \defgroup Conditions_limites
Classe relative aux conditions limites
*/

using namespace LMT;
using namespace Codegen;
//******************************************
// Condition aux limites
//*******************************************
/** \ingroup Conditions_limites
\brief Classe definissant les conditions aux limites
*/
struct Boundary {
    static void buildTimeSymbols(const DataUser &data_user);                        ///< Construit la liste des symboles temporels
    static void updateTimeValues(const Process &process,const DataUser &data_user); ///< Actualise les valeurs temporelles
    static void buildSpaceSymbols();                                                ///< Construit la liste des symboles spatiaux
    static void updateSpaceValues(const Vec<TYPEREEL,DIM> &space_point);            ///< Actualise les valeurs spatiales aux coordonnees du point "space_point"
    static void buildBoundarySymbols(const DataUser &data_user);                    ///< Construit les listes de symboles permettant d'evaluer les CL
    static void free();
    
    /// Donnees communes aux CL
    static std::vector<Ex> time_symbols;   ///< Vecteur contenant les symboles temporels
    static Ex::MapExNum time_values;        ///< Map< Symbol, Valeur > des grandeurs symboliques (autres que spatiales) du calcul
    static std::vector<Ex> space_symbols;  ///< Vecteur contenant les symboles spatiaux
    static Ex::MapExNum space_values;       ///< Map< Symbol, Valeur > des grandeurs symboliques spatiales du calcul
    
    Boundary() {ft.set(0.0);}
    ~Boundary() {free();}
    
    /** Evalue la valeur de la CL en chaque noeud de "nodeeq" au step numero "i_step" et stocke le resultat dans V
     * ATTENTION!!! Cette fonction ne reactualise pas les valeurs dependant du temps!
     */
    void evaluate(unsigned i_step, Vec<TYPEREEL> &V, Vec<Vec<TYPEREEL,DIM> > &nodeeq,const Vec<TYPEREEL> &neqs = Vec<TYPEREEL>());
    
    int id;                                     ///<  id de la condition limite, la même que dans data_user
    Sc2String comp;                             ///< type de condition aux limites
    Vec<Vec<TYPEREEL,DIM>, 2> box, box1;        ///< boite incluant la CL dans laquelle sont cherchees les interfaces, box1 la boite de la 2eme interface pour les périodiques
    Vec<Vec<Sc2String,DIM> > fcts_spatiales;    ///< fonctions chargement spatial : 2 ou 3 fcts analytiques (x, y, z) separees par des ; (premier vecteur égal au nombre de step)
    Vec<Sc2String > fcts_temporelles;           ///< fonctions chargement temporel : 1 fct de "t" par step pour chaque coordonnée      
    Vec<Vec<TYPEREEL,2> > intervalles_temps;    ///< intervalles des differentes fonctions temporelles
    Vec<TYPEREEL,DIM> ft;                       ///< valeur de la fonction temporelle pour le piquet de temps considere pour chaque coordonnee
    Vec<unsigned> sst_num;                      ///< numero des sous-structures où chercher la CL
};

#endif //BOUNDARY_H
