// Author: D. Violeau , 24/03/2006

//librairies Hugo
#include "containers/mat.h"
#include "containers/vecpointedvalues.h"
#include "mesh/mesh.h"
#include <fstream>
#include <map>
#include "containers/vec_mt.h"

//fichiers de definition des variables
#include "definition_CL.h"
#include "definition_PARAM.h"
#include "definition_PARAM_STRUCTURE.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"

#include "assignation_mpi.h"

#include "create_geometry_mesh.h"

using namespace LMT;
using namespace std;

/** \defgroup Maillage_geometrie Création des sous-structures et interfaces
\brief Génération des maillages et géométrie des sous-structures, interfaces.
 
La fonction principale est multiscale_geometry_mesh() écrite en template. 
La fonction fake_multiscale_geometry_mesh() permet de forcer le compilateur à compiler avec les types voulus la fonction multiscale_geometry_mesh()

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


void fake_multiscale_geometry_mesh() {
    XmlNode n;
    Param process;

    Vec<Sst<DIM,TYPEREEL> > S3;
    Vec<VecPointedValues<Sst<DIM,TYPEREEL> > > SubS3, Stot3;
    Vec<Interface<DIM,TYPEREEL> > Inter3;
    Vec<VecPointedValues<Interface<DIM,TYPEREEL> > > SubI3;
    Vec<Boundary<DIM,TYPEREEL> > CL3;
    
    multiscale_geometry_mesh(n,S3, Inter3, process, CL3,Stot3,SubS3,SubI3) ;
}


