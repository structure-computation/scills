using namespace LMT;
using namespace std;

#include "assign_material.h"

//********************************************
// Affectation de materiaux aux SST
//******************************************** 
/**\ingroup Materiaux
 \brief Structure contenant un operateur() (utilisable avec apply) permettant l'affectation des valeurs des proprietes materiaux au maillage associé à la sous-structure.
 
 Chaque sous-structure comporte un champ typmat (Sst::typmat) dont la valeur correspond à une composante du vecteur des propriétés matériaux SstCaract. 
 Selon le nom correspondant extrait du vecteur de propriétés matériaux, on déclare une formulation du type donné. 
 Si certaines formulations ne veulent pas être inclus dans la base des formulations disponibles, il suffit de commenter chaque bloc if et d'enlever la formulation correspondante dans le fichier SConstruct.  
 */
struct assignation_material_to_SST{
   template<class SST, class TV3> void operator()(SST &S,TV3 &matprop, Param &process,DataUser &data_user) const{
//       S.mesh.load();
      assign_material(S,matprop,process,data_user); //fichier genere par __init__.py dans LMT
//       S.mesh.unload();
   }
};
