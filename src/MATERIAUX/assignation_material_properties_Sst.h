#ifndef ASSIGNATION_MATERIAL_PROPERTIES_SST_H
#define ASSIGNATION_MATERIAL_PROPERTIES_SST_H


#include "MATERIAUX_declarations.h"


/** \ingroup Materiaux
\brief Lecture des proprietes materiau par sst.

Le bloc <materials>  </materials> est lu dans le fichier xml.

On détermine ensuite le nombre de blocs <coefficients /> pour créer le vecteur de propriétés matériau (classe SstCarac).

L'"identificateur" permet de positionner le matériau lu dans le vecteur de propriétés. Le numéro affecté aux sous-structures correspond alors à cet identificateur.

Selon la dimension le champ "resolution" est lu ou non. Celui ci peut prendre les valeurs
\code resolution="deformation_plane" ou resolution="contrainte_plane" \endcode

Ensuite le champ "type" renseigne sur le comportement de la sous-structure. Selon la valeur obtenue, les différentes propriétés sont extraites du fichier xml (les champs nécessaires sont consultables dans les fichiers de formulations correspondants).
- \code type="isotrope" \endcode
- \code type="orthotrope" \endcode
- \code type="viscoelas" \endcode
- \code type="orthotrope_damage" \endcode

On aura donc les possibilités suivantes : 
-   On peut choisir un matériau isotrope  :
\code 
<coefficients type="isotrope" resolution="contrainte_plane" identificateur="0" name="fibre"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6"/>
\endcode
-   ou bien un matériau orthotrope
\code
<coefficients type="orthotrope" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0"/>
\endcode
-   ou encore un matériau orthotrope endommageable pour lequel on ajoute 
\code
<coefficients type="orthotrope_damage" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0" Yo="0.0961" Yop="0.0961" Ycp="10.8" Yc="10.8" b="2.5" Ysp="0.70"/>
\endcode
-   ou encore un matériau viscoelastique (orthotrope ou isotrope)
\code
<coefficients type="isotrope" resolution="contrainte_plane" identificateur="0"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6" viscosite="0.1" />
\endcode
*/
void read_matprop(Vec<SstCarac> &matprops, Param &process, DataUser &data_user , BasicVec<BasicVec<TYPEREEL> > &mat_prop_temp);


/** \ingroup Materiaux
\brief Structure contenant un operateur() (utilisable avec apply) permettant l'affectation des valeurs des proprietes materiaux au maillage associé à la sous-structure.

Chaque sous-structure comporte un champ typmat (Sst::typmat) dont la valeur correspond à une composante du vecteur des propriétés matériaux SstCaract. 
Selon le nom correspondant extrait du vecteur de propriétés matériaux, on déclare une formulation du type donné. 
Si certaines formulations ne veulent pas être inclus dans la base des formulations disponibles, il suffit de commenter chaque bloc if et d'enlever la formulation correspondante dans le fichier SConstruct.  
*/
struct assignation_material_to_SST{
    void operator()(Sst &S,Vec<SstCarac> &matprops) const;
};


#endif //ASSIGNATION_MATERIAL_PROPERTIES_SST_H