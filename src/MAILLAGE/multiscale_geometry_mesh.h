#ifndef MULTISCALE_GEOMETRY_MESH_H
#define MULTISCALE_GEOMETRY_MESH_H

#include "../all_declarations.h"

/** \defgroup Maillage_geometrie Création des sous-structures et interfaces
 \ brief Génération des maillages et géométrie des sous-structures, inte*rfaces.
 
 La fonction principale est multiscale_geometry_mesh() : onction principale de création des maillages et géométries (voir ci-dessous). 
 
 L'ordre des opérations est le suivant :
 - create_SST_INTER() : lecture et création des sst et interfaces
 - create_SST_typmat() : Création du vecteur S de sst et assignation d'un matériau à partir du fichier de qualification
 - create_maillage_SST() : Lecture et assignation des maillages des ssts
 - create_perfect_interfaces() : Création des interfaces entre ssts
 - create_gap_interfaces() : Création des interfaces avec jeu physique
 - create_interfaces_CL() : Création des interfaces comprises dans les CL
 - calcul_taille_probleme() : Calcul de la taille du problème
 - calculate_measure_G_neq_INTER : calcul de la mesure de G et des normales équivalentes pour l'interface
 - barycenter_constant_rho() : Calcul du centre de gravité G d'un maillage
 - measure() : Calcul de la mesure à partir d'un maillage
 - calculate_normales : Calcul des normales de chaque coté de l'interface
 - selon le choix de sous-integration :
 - surdiscretise_SST et decoup_INTER: sur-discrétisation h des maillages des ssts, création des maillages de bord des sst
 - p_surdiscretise_SST et p_decoup_INTER : sur-discretisation p des maillages des sst, création des maillages de bord des sst
 - copy_INTER : création des maillages de bord des sst
 - calculate_measure_G_SST : calcul de la mesure et cdg des sst
 
 */



/** \ingroup Maillage_geometrie
 \ brief Fonction principale de création des maillages et géométries    *
 
 - Dans un premier temps on crée les sous-structures et interfaces à partir des informations contenues dans la class GeneralParameters (inclus dans Process) et des conditions aux limites CL. On obtient ainsi un vecteur de sous-structures \ref Sous_structures et d'interfaces \ref Interfaces.
 La creation des sous-structures et interfaces est effectuée de la facon indiquée dans la fonction create_SST_INTER() .
 
 On crée ensuite les données géométriques associées aux interfaces : calculate_measure_G_neq_INTER.
 
 - Dans un second temps, on modifie le maillage des sous-structures pour prendre en compte ou non (d'après le paramètre Process::sousint) la sous-integration (cf. surdiscretise_SST ).
 
 - Enfin on détermine certaines caractéristiques géométriques des Sst (Sst::G, Sst::measure, Sst::Edge::G).
 
 Rq : La réorganisation des noeuds des Sous-structures n'est pas nécessaire étant donné qu'on utilise une renumérotation dans le solveur sparse (symamd).
 
 */
void multiscale_geometry_mesh(DataUser                          &data_user, 
                              GeometryUser                      &geometry_user,
                              Vec<Sst>                          &S,
                              Vec<Interface>                    &Inter, 
                              Process                           &process, 
                              Vec<Boundary>                     &CL,
                              Vec<VecPointedValues<Sst> >       &Stot,
                              Vec<VecPointedValues<Sst> >       &SubS,
                              Vec<VecPointedValues<Interface> > &SubI);

#endif //MULTISCALE_GEOMETRY_MESH_H