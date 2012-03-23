#ifndef ASSIGNATION_MATERIAL_PROPERTIES_INTERFACE_H
#define ASSIGNATION_MATERIAL_PROPERTIES_INTERFACE_H


#include "MATERIAUX_declarations.h"


/** \ingroup Materiaux
\brief Lecture des proprietes materiau d'interface 
 
On recherche le bloc <proprietes_interfaces>  </ proprietes_interfaces> puis on détermine le nombre de propriétés différentes en considérant le nombre de champs <coefficients />. On crée ainsi un vecteur de propriétés d'interface.
 
On a la possibilité d'entrer différentes propriétés particulières : 
- interface de type contact renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_sst" coeffrottement="0.3" num_sst="0 1" name="cube_cube"/>
\endcode
- interface de type contact renseignée par une boite déterminée par ses deux points extrèmes (classés par ordre croissant).
\code 
<coefficients type="contact_box" coeffrottement="0.3" box="1 0 0 1 10 10" name="cube_cube"/>
\endcode
Le champ box contient donc le point inférieur gauche puis le point supérieur droit (4 composantes en 2d 6 en 3d). Toutes les interfaces dont la boite est incluse dans cette boite sont repérées. Attention pour les interfaces courbes (la boite n'est plus un plan ou une droite).
- interface de contact avec jeu renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_jeu_sst" coeffrottement="0.3" num_sst="0 1" jeu="x*x+y+1" name="cube_cube"/>
\endcode
Le jeu est une fonction quelconque des coordonnées des interfaces concernées dans le repère x, y, z. S'il est rentré comme précédement, il sera considéré comme étant normal aux surfaces. On peut également le rentrer complet <tt>jeu="x*x;0;0"</tt> ce qui construira un champs par point fonction de l'espace.
- interface de contact avec jeu renseignée par une boite déterminée par ses deux points extrêmes (classés par ordre croissant).
\code 
<coefficients type="contact_jeu_box" coeffrottement="0.3" box="0 0 0 0 0 0" jeu="x*x+y+1" name="fibre-matrice"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de type cohésive renseignée par une boite.
\code 
<coefficients type="cohesive" coeffrottement="0.3" kn="0.12" kt="0.12" knc="0.12" gamma="0" alpha="0" Yc="0" Yo="0" n="2.5" name="fibre-matrice"/>
\endcode
- interface de type discrète
\code
<coefficients type="discrete" coeffrottement="0.3" Gcrit="0.12"  name="fibre-matrice"/>
\endcode
Lorsque le taux de restitution critique est inférieur à une valeur critique cette interface est parfaite sinon elle est de type contact avec frottement. Le nom est ici important puisque c'est à partir de celui-ci que les propriétés matériaux sont assignées (on recherche de quel type sont les sous-structures adjacentes)
- interface de type jeu imposé renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="jeu_impose_sst" num_sst="0 1" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de contact avec jeu renseignée par une boite déterminée par ses deux points extrêmes (classés par ordre croissant).
\code 
<coefficients type="jeu_impose_box" box="0 0 0 0 0 0" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
*/
void read_propinter(Vec<InterCarac> &propinter,const DataUser &data_user, BasicVec<BasicVec<TYPEREEL> > &link_prop_temp);

/** \ingroup Materiaux
 \ *brief Assignation des propriétés aux interfaces
 
 On alloue tout d'abord la mémoire pour les paramètres matériaux (cf. PARAM_COMP_INTER) de toutes les interfaces.
 
 On boucle ensuite sur les propriétés d'interfaces (Matprop) et selon le type (box ou sst) on sélectionne les interfaces correspondantes.
 
 Enfin pour ces interfaces sélectionnées seulement, on modifie le champ Interface::comp, on assigne ensuite les paramètres aux interfaces.
 
 Pour le jeu, en utilisant codegen et les noeuds équivalents de l'interface, on crée le champ jeu qui sera par la suite utilisé dans l'étape locale. Si l'utilisateur renseigne une seule valeur, on considère que c'est un jeu selon la normale, sinon il faut renseigner les valeurs selon x, y (et z) du jeu. 
 */
void modif_inter(Vec<Interface> &Inter, Vec<InterCarac> &propinter, Param &process,const DataUser &data_user);


#endif //ASSIGNATION_MATERIAL_PROPERTIES_INTERFACE_H