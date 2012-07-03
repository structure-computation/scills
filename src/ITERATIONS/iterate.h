#ifndef ITERATE_H
#define ITERATE_H

// fonction utilisees pour la procedure iterative
#include "../DEFINITIONS/Process.h"
#include "../DEFINITIONS/Structure.h"
#include "../DEFINITIONS/MacroProblem.h"

using namespace LMT;


/** \defgroup Incrementale Strat�gie Incr�mentale
\ingroup Strategie_iterative
\brief Proc�dure principale de la strat�gie incr�mentale en quasistatique.
 
 
 Pour chaque piquet de temps, on effectue une r�solution it�rative. On assigne une valeur diff�rente au coefficient LatinData::mu pour la relaxation. A la premi�re it�ration on lui assigne 1.0 sinon on lui affecte la valeur de LatinData::facteur_relaxation (donn� par l'utilisateur).
 
 On effectue la boucle latin tant que le crit�re d'erreur n'est pas atteint (LatinData::critere_erreur) ou pour un nombre maximal d'it�rations (LatinData::nbitermax). Celle ci est d�crite dans les modules suivants :
- \ref etape_lineaire pour le piquet de temps 1
- \ref relaxation_quantites pour le piquet de temps 1
- \ref etape_locale pour le piquet de temps 1
- \ref calcul_erreur_incr
 
 Le param�tre LatinData::list_error permet de lister l'erreur latin au cours des it�rations.
*/
void iterate_incr(Process &process, PointedSubstructures &S, VecInterfaces &Inter,PointedInterfaces &SubI, MacroProblem &Global,DataUser &data_user);

#endif // ITERATE_H
