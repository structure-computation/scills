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

/** \defgroup Maillage_geometrie Cr�ation des sous-structures et interfaces
\brief G�n�ration des maillages et g�om�trie des sous-structures, interfaces.
 
La fonction principale est multiscale_geometry_mesh() �crite en template. 
La fonction fake_multiscale_geometry_mesh() permet de forcer le compilateur � compiler avec les types voulus la fonction multiscale_geometry_mesh()

L'ordre des op�rations est le suivant :
- create_SST_INTER() : lecture et cr�ation des sst et interfaces
   - create_SST_typmat() : Cr�ation du vecteur S de sst et assignation d'un mat�riau � partir du fichier de qualification
   - create_maillage_SST() : Lecture et assignation des maillages des ssts
   - create_perfect_interfaces() : Cr�ation des interfaces entre ssts
   - create_gap_interfaces() : Cr�ation des interfaces avec jeu physique
   - create_interfaces_CL() : Cr�ation des interfaces comprises dans les CL
   - calcul_taille_probleme() : Calcul de la taille du probl�me
- calculate_measure_G_neq_INTER : calcul de la mesure de G et des normales �quivalentes pour l'interface
   - barycenter_constant_rho() : Calcul du centre de gravit� G d'un maillage
   - measure() : Calcul de la mesure � partir d'un maillage
   - calculate_normales : Calcul des normales de chaque cot� de l'interface
- selon le choix de sous-integration :
   - surdiscretise_SST et decoup_INTER: sur-discr�tisation h des maillages des ssts, cr�ation des maillages de bord des sst
   - p_surdiscretise_SST et p_decoup_INTER : sur-discretisation p des maillages des sst, cr�ation des maillages de bord des sst
   - copy_INTER : cr�ation des maillages de bord des sst
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


