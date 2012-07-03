#ifndef ITERATE_H
#define ITERATE_H

// fonction utilisees pour la procedure iterative
#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Structure.h"
#include "../DEFINITIONS/MacroProblem.h"

using namespace LMT;


/** \defgroup Incrementale Stratégie Incrémentale
\ingroup Strategie_iterative
\brief Procédure principale de la stratégie incrémentale en quasistatique.
 
 
 Pour chaque piquet de temps, on effectue une résolution itérative. On assigne une valeur différente au coefficient LatinData::mu pour la relaxation. A la première itération on lui assigne 1.0 sinon on lui affecte la valeur de LatinData::facteur_relaxation (donné par l'utilisateur).
 
 On effectue la boucle latin tant que le critère d'erreur n'est pas atteint (LatinData::critere_erreur) ou pour un nombre maximal d'itérations (LatinData::nbitermax). Celle ci est décrite dans les modules suivants :
- \ref etape_lineaire pour le piquet de temps 1
- \ref relaxation_quantites pour le piquet de temps 1
- \ref etape_locale pour le piquet de temps 1
- \ref calcul_erreur_incr
 
 Le paramètre LatinData::list_error permet de lister l'erreur latin au cours des itérations.
*/
void iterate_incr(Process &process, PointedSubstructures &S, VecInterfaces &Inter,PointedInterfaces &SubI, MacroProblem &Global,DataUser &data_user);

#endif // ITERATE_H
