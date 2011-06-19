 #ifdef METIL_COMP_DIRECTIVE
// #pragma src_file LMT/include/util/solveLDL.cpp
// #pragma src_file LMT/include/util/solveLDL_semimorse.cpp
// #pragma src_file LMT/include/util/symamd.cpp
// #pragma src_file LMT/include/util/symrcm.cpp
// #pragma cpp_path /usr/include/suitesparse/
// #pragma inc_path /usr/include/suitesparse/
// #pragma lib_name cholmod
// 
// #pragma inc_path /home/ubuntu/driver_toolkit/NVIDIA_GPU_Computing_SDK/C/common/lib /home/ubuntu/driver_toolkit/NVIDIA_GPU_Computing_SDK/C/lib/
// #pragma lib_name lam 
// #pragma lib_name ldl 
// #pragma gpu_flag --maxrregcount 62 -arch=sm_21 

#pragma src_file multiscale_geometry_mesh.cpp
#pragma src_file assignation_materials_property.cpp
#pragma src_file multiscale_operateurs.cpp
#pragma src_file formulation_2_double_elasticity_isotropy_stat_Qstat.cpp
#pragma src_file formulation_3_double_elasticity_isotropy_stat_Qstat.cpp
#pragma src_file formulation_2_double_elasticity_orthotropy_stat_Qstat.cpp
#pragma src_file formulation_3_double_elasticity_orthotropy_stat_Qstat.cpp
#pragma src_file iterate_stat_Qstat.cpp
#pragma src_file affichage.cpp

#endif

//librairie Hugo
#include "containers/mat.h"

// fichiers de definition des classes
#include "definition_PARAM.h"
#include "definition_GLOB.h"
#include "definition_SST_time.h"
#include "definition_INTER_time.h"
#include "definition_CL.h"

//allocation memoire
#include "allocation_memoire_PARAM.h"

//lecture des parametres
#include "multiscale.h"

//assignation des différents paramètres MPI
#include "assignation_mpi.h"


// #include "containers/evaluate_nb_cycles.h"
#include "SCtime.h"

//biblioteque venant de SC_create_2
#include <Metil/BasicVec.h>
#include <Metil/StructCompactor.h>
#include <Metil/CudaMetil.h>
#include <Metil/Hdf.h>
using namespace Metil;

#include "GeometryUser.h"
#include "DataUser.h"

#ifndef INFO_TIME
#define INFO_TIME
#endif 

using namespace LMT;
using namespace std;
Crout crout;
#ifdef PRINT_ALLOC
namespace LMT {
    std::map<std::string,long long> total_allocated;
};
#endif

//*********************************************************************
// procedure multiechelle quasistatique ou statique (2 ou 3d)
//*********************************************************************
int main(int argc,char **argv) {

    /// lecture du fichier de donnees
    if ( argc!=3 and argc !=4 ) {
        std::cerr << "usage : nom_executable + dimension +  nom_complet_du_fichier_de_donnees. (+ mpi)" << std::endl;
        std::cerr << "ex : ./SC_multi_2 num_model num_calcul  (+ mpi)" << std::endl;
        return 1;
    }
    try {
        Param process;
        allocation_memoire_PARAM(process);
        /// on fait ce qu il y a faire quand on est en MPI
        if (argc == 4 and strcmp(argv[3], "mpi")==0 ) {
            definition_mpi_param(process,argc,argv);
            if (process.rank == 0) std::cout << "Calcul MPI" << std::endl;
        } else {
            process.rank=0;
            process.size=1;
        }
        crout.open(process.rank);
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    TicTac tic1;
    if (process.rank==0) {tic1.init();tic1.start();}
#endif
        
        if (process.rank == 0 ) std::cout << "*********************************************" << std::endl;
        if (process.rank == 0 ) std::cout << "* Lecture des données venant de SC-create_2 *" << std::endl;
        if (process.rank == 0 ) std::cout << "*********************************************" << std::endl;
        if (process.rank == 0 ) std::cout << std::endl;
    
        std::string id_model = argv[ 1 ];
        std::string id_calcul = argv[ 2 ];
        
        // ******************************************************************************************************************
        //lecture des données utilisateur (fichier de calcul .json) et compilation à la volée 
        DataUser data_user(id_model, id_calcul);
        data_user.read_json_calcul();

        std::cout << "lecture de la geometrie" << std::endl;
        // ******************************************************************************************************************
        //lecture de la geometrie--------------------------------------------------
        GeometryUser geometry_user(id_model, id_calcul);
        geometry_user.read_hdf5(false,true);                       // true si on lit les info micro, true si on lit toutes les infos
        std::cout << "fin lecture de la geometrie" << std::endl;
        geometry_user.split_group_edges_within_geometry(data_user);
        std::cout << "fin lecture de la geometrie" << std::endl;
    
        if (process.rank == 0 ) std::cout << "*****************************" << std::endl;
        if (process.rank == 0 ) std::cout << "* DECOMPOSITION DE DOMAINES *" << std::endl;
        if (process.rank == 0 ) std::cout << "*****************************" << std::endl;
        if (process.rank == 0 ) std::cout << std::endl;


        //Creation du dossier tmp
        system("mkdir -p tmp");

        if (process.rank == 0 ) std::cout << "****************************" << std::endl;
        if (process.rank == 0 ) std::cout << " Lecture des donnees " << std::endl;
        if (process.rank == 0 ) std::cout << "****************************" << std::endl;

        /// lecture des donnees de calcul
        // donnees associees au calcul
        process.dim=DIM;
        process.affichage->name_data= argv[2];
        if (process.rank == 0 ) std::cout << "\tFichier lu : " << process.affichage->name_data << std::endl;

//         Vec<Splitted<Sst<DIM,TYPEREEL> , 16> > S;
        Vec<Sst<DIM,TYPEREEL> > S;
        Vec<Interface<DIM,TYPEREEL> > Inter;
        Vec<Boundary<DIM,TYPEREEL> > CL;
        Glob<DIM,TYPEREEL> Global;


        
        S.resize(geometry_user.nb_group_elements);
        multiscale(data_user, geometry_user, S, Inter, process,  CL, Global);

        
#ifdef INFO_TIME
    if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
    if (process.rank==0) std::cout << "Duree complete du programme : " ;
    if (process.rank==0) tic1.stop();
    if (process.rank==0) std::cout << std::endl;
#endif
        if (process.size>1) MPI_Barrier(MPI_COMM_WORLD);
        desallocation_memoire_PARAM(process);

	if (process.size > 1 ) MPI_Finalize();

    } catch ( const IoException &e ) {
        std::cout << e.what() << std::endl;
    }
    return 0;
}




/** \mainpage Fichier general pour le programme multiechelle.
 
Ce programme est constitué de différentes sections chacune correspondant à un fichier cpp séparé :
 
- \ref lecture "Lecture des données"
 
- \ref Maillage_geometrie "Création des sous-structures et interfaces"
 
- \ref Materiaux "Assignation des matériaux aux sous-structures et des propriétés aux interfaces"
 
- \ref Operateurs "Création des opérateurs"
 
- \ref Strategie_iterative "Stratégie itérative"
 
- \ref Exemples "Exemples de fichiers de données xml"
*/


/** \page lecture Lecture des donnees.
 
Dans un premier temps le fichier xml entré comme argument est lu et 
les différentes quantités sont affectées aux instances des classes correspondantes. 
 
L'ordre dans lequel sont rentrées les différents champs xml n'a pas d'importance.
 
Les champs du fichier xml doivent être inclus dans  :
 
\code  
<?xml version="1.0"?>
<root>
 
</root>
\endcode
 
\section Lecture_calcul Lecture des données liées au calcul.
La lecture du xml permet d'assigner les différents champs de l'instance process de la class Param.
Les champs xml suivants sont lus (cliquer sur les liens pour voir le détail des sous-champs à remplir):
- paramètres assignés à la classe MULTI et Param \ref parametres_xml    
\code <parametres /> \endcode               
- paramètres temporels assignés à la classe TEMPS. \ref parametres_temporels_xml
\code <parametres_temporels /> \endcode  
- parametres d'affichage (classe AFFICHAGE). \ref parametres_affichage_xml
\code <parametres_affichage /> \endcode  
- paramètres pour les directions de recherche (classe MULTI). \ref parametres_recherche_xml 
\code <parametres_recherche /> \endcode  
 
\section Lecture_structure Lecture des données associées à la géométrie et au maillage 
On lit ensuite la ligne suivante correspondant aux maillages et nombre de maillages à lire que l'on assigne à l'instance de la classe STRUCTURE contenue dans Param :\ref mesh
\code 
    <mesh />
\endcode
 
\section Lecture_CL Lecture des conditions aux limites.
On termine par la lecture des différentes lignes de conditions aux limites (1 bloc CL par condition ce qui correspond à une entrée dans le vecteur de conditions aux limites classe Boundary) : \ref CL
\code 
   <CL>
   
   </CL>
\endcode
 
 
 
 
\section Details Détails des champs xml
 
\subsection parametres_xml - Détails parametres
Liste exhaustive des champs disponibles :
\code <parametres sous_integration="1" type_sous_integration="p" nbmacro="9" multiechelle="1" blocage_modes_rigides="0" mvts_bloques="Ty Tx Rz" nbitermax="800" facteur_relaxation="0.8" critere_erreur="1e-5" type_erreur="ddr" nb_threads="1" calcul_micro="non" deltaT="0"/>   \endcode
- sous-integration (obligatoire) : 1 pour utiliser la sous-integration 0 sinon.
- type_sous_integration (obligatoire) : dans le cas ou sous-integration=="oui" on lit ce paramètre qui prend la valeur "h" pour un redécoupage des éléments (ne marche que pour les triangles quad et tetra) ou la valeur "p" pour une augmentation du degré d'interpolation de l'élément (passage de P1 à P2 par ajout de noeuds partout)
- nbmacro (obligatoire) : nombre de fonctions macro spatiales pour toutes les interfaces (9 en 3D et 4 en 2D).
- multiechelle (obligatoire) : booleen = 1 pour effectuer un calcul multiechelle, 0 pour un calcul monoechelle
- blocage_modes_rigides (obligatoire) : booléen =1 pour dire de bloquer explicitement les modes de corps rigides.
- mvts_bloques (obligatoire) : vecteur contenant les valeurs Tx Ty Tz Rx Ry Rz (séparés par des espaces). A renseigner uniquement si blocage_modes_rigides == 1. Par défaut lorsqu'il n'y a que des interfaces de type effort imposé sur le bord, le programme détermine automatiquement les mouvements macro à bloquer. Par contre s'il ne faut bloquer qu'une translation (exemple lorsqu'on utilise une condition de symétrie), il est nécessaire d'aider le programme en indiquant la translation à bloquer.
- nbitermax (obligatoire) : nombre d'itérations maximal pour la LATIN. Si le critère est atteint avant on sort de la boucle sur les itérations.
- facteur_relaxation (obligatoire) : 0.8 est souvent une bonne valeur
- critere_erreur (obligatoire) : critère à atteindre pour arréter l'algorithme latin.
- type d'erreur (obligatoire) : choix possible entre "ddr" et "energie" (cf. \ref calcul_erreur ).
- nb_threads (obligatoire) : nombre de processus à lancer sur une même machine
- deltaT (obligatoire) : valeur du champ de température pour un calcul thermoélastique
 
\subsection parametres_temporels_xml - Détails parametres_temporels
Liste exhaustive des champs disponibles :
\code <parametres_temporels type_de_calcul="Qstat" nbpastemps="24" pasdetemps="0.05" /> \endcode
- type_de_calcul (obligatoire) : choix entre "Qstat" pour la quasistatique et "stat" pour la statique.
- nbpastemps (non nécessaire en statique) : nombre d'intervalle d'étude (ou piquet de temps)
- pasdetemps (non nécessaire en statique) : valeur du pas de temps
 
\subsection parametres_affichage_xml - Détails parametres_affichage
Liste exhaustive des champs disponibles :
\code   <parametres_affichage interactivite="1" affich_mesh="0" affich_resultat="1" type_affichage="Sinterieur" save_or_display="display" display_fields="dep qtrans num sigma epsilon ener" list_error="1" display_error="1" repertoire_save="./tmp/3d_2_cubes/" name_data="data1" affich_depl_pt="1" coor_point="100 100 100"/> \endcode 
- interactivite (obligatoire) : booléen = 1 permettant à la fin d'un calcul de post-traiter les résultats (cf. Interactivité), 0 sinon et par défaut on utilise alors les autres champs pour afficher le résultat ou non.
- affich_mesh (obligatoire) : (booléen =1 ) une fois les maillages créés et surdiscrétisés pour les ssts et les interfaces, il est possible de visualiser les maillages avec éclaté des sous-structures ou interfaces, type de matériau ou d'interface et numéro.
- type_affichage (obligatoire) : choix parmi "Sinterieur", "Sbord" et "Inter" pour afficher le maillage des sous-structures ou de bord des ssts ou des interfaces.
- save_or_display (obligatoire) : choix parmi "save" ou "display" : save pour ne sortir que le fichier paraview0.vtu ou display pour le créer et lancer paraview.
- affich_resultat (obligatoire) : (booléen = 1) affichage du résultat après critère d'erreur atteint. On crée alors un fichier pvd permettant de visualiser la déformée ou autre au cours des pas de temps.
- display_fields (obligatoire) : champs séparés par des espaces permettant de sélectionner les champs à afficher. Pour la liste des noms disponibles, se reporter aux fichiers de formulations dans le répertoire FORMULATIONS. 
- repertoire_save (obligatoire) : répertoire de sauvegarde des résultats de calcul.
- name_data (obligatoire) : nom générique pour chaque pas de calcul (ex : pour le pas 2 on aura repertoire_save/data1_2.vtu). Le fichier pvd créé sera lui aussi dans le répertoire sous le nom : repertoire_save/data1.pvd
- list_error (obligatoire) : booléen = 1 pour faire défiler l'erreur au cours des itérations
- display_error (obligatoire) : booléen = 1 pour afficher la courbe d'erreur en fin de calcul
 
\subsection parametres_recherche_xml - Détails parametres_recherche
Liste exhaustive des champs disponibles :
\code  <direction_recherche ktype="scalaire_auto_CL" kfact="1" formule="E/L" copydirection="1" /> \endcode
- ktype (obligatoire) : type de direction de recherche : les noms disponibles sont détaillés dans LATIN::ktype, le mieux est d'utiliser "scalaire_auto_CL" avec kfact=1
- kfact (obligatoire) : facteur multiplicatif des directions de recherche ou directement la direction de recherche
- copydirection (obligatoire) : booleen =1 : on effectue une moyenne des directions obtenues que l'on assigne aux deux cotés de l'interface
 
\subsection mesh - Détails mesh
Liste exhaustive des champs disponibles :
\code  <mesh repertoire_des_maillages="./DATAFILES/test3D/3d_isotrope_cubes/" nom_fichier_qualification_materiaux="Vqualif.py" nom_des_maillages="V" nb_de_maillages="2" extension=".avs" jeu_physique="0" /> \endcode
- repertoire_des_maillages (obligatoire) : contient les maillages des sous-structures à lire
- nom_des_maillages (obligatoire) : nom générique pour les maillages (ex : V0.avs, V1.avs ...). Le premier numéro est nécessairement 0.
- extension (obligatoire) : permet de lire des fichiers différents selon l'extension
- nb_de_maillages (obligatoire) : correspond au nombre de sous-structures
- jeu_physique (obligatoire) : booléen = 1 pour lire les maillages d'interfaces disjointes
- nom_maillages_jeu (uniquement si jeu_physique==1) : nom générique des maillages d'interfaces donnés. Les noms doivent être générés de la facon suivantes : nom_maillages_jeu0-0-2.avs -> interface extraite de la sst numero 0 et ayant comme voisins les sst 0 et 2.
- inter_jeu (uniquement si jeu_physique==1) : vecteur de paire de sous-structures séparés par des ; (ex : "0 2;1 3")
 
 \subsection CL - Détails CL
 
 A l'intérieur d'un bloc plusieurs champs sont disponibles selon le type de CL
\code     
<CL>
      <parametres box="100 0 100    100 100 200" unit="mm" comp="effort" />
      <fct_spatiale fonction="1;0;0" />
      <fct_temporelle intervalle="0. 0.1"  fonction="0" />
      <fct_temporelle intervalle="0.1 1.1"  fonction="-50*sin(3.14159265358979*(t-0.1))" />
   </CL>
 \endcode 
   - le champ box contient une boite définie par ses deux extrémités (2*2 composantes en 2d et 3*2 en 3d) permettant d'extraire le maillage inclue dans cette boite
   - comp renseigne sur le type de comportement (effort, depl, sym, depl_normal)
   - le champ fct_spatiale n'est utilisé que pour les interfaces de type effort ou deplacement. A savoir, on rentre les valeurs des fonctions imposées selon x, y (et z en 3d) séparées par des ; (ne pas mettre d'espace). Les fonctions doivent être des fonctions de x, y (et z en 3d) selon la zone considérée. Ex : surface (0 0  0 1 ) ->  fonction de y
   Pour une condition de type deplacement normal imposé (depl_normal) seule une fonction est nécessaire, on doit aussi ajouter le champ normale (x,y ou z). Pour une condition de symétrie ce champ peut être effacé.
   \code 
   <fct_spatiale fonction="(0*y/2+0)" normale="x"/>
   \endcode
   - On renseigne autant de blocs fct_temporelle que d'intervalles pour lesquels la fonction temporelle est différente. La valeur de la fonction est évaluée pour chaque pas de temps et multipliée à la fonction spatiale.
   
   Conseil : Il est préférable de prendre une fct_spatiale unitaire pour que le produit fct_spatiale*fct_temporelle soit simple à interpréter.
 
 */



/** \page Exemples Exemples de fichiers de données xml
 
Voici une liste d'exemples de fichiers de données xml compris dans le répertoire DATAFILES
 
Le premier exemple est une structure 2D chargée uniquement par des efforts ou déplacements (la variation de température est nulle). On présente ici deux possibilités : 
 - contrainte plane avec chargement en effort :  \ref data_stat_isotrope_CP_effort.xml
 - déformation plane avec chargement normal (déformation longitudinale) : \ref data_stat_isotrope_DP_depl.xml
 
Le deuxième exemple est une structure 3D orthotrope (2 matériaux) avec chargement en température uniquement : \ref data_stat_ortho_therm.xml
 
Le troisième exemple est une structure 3D isotrope en quasistatique avec évolution de l'effort imposé : \ref data_Qstat_iso.xml
 
Le quatrième exemple est une structure 3D isotrope avec des interfaces de type contact avec frottement et sans jeu en quasitatique avec évolution du chargement : \ref data_Qstat_iso_contact.xml
 
Le dernier exemple est une structure 3D isotrope avec des interfaces de type jeu (contact avec frottement) en quasistatique avec évolution du chargement :  \ref data_Qstat_iso_jeu.xml
 
*/

/** \example data_stat_isotrope_CP_effort.xml
\brief Structure 2D chargée uniquement par des efforts ou déplacements (la variation de température est nulle) : contrainte plane avec chargement en effort
*/

/** \example data_stat_isotrope_DP_depl.xml
\brief Structure 2D chargée uniquement par des déplacements normaux (la variation de température est nulle) : déformation plane
*/

/** \example data_stat_ortho_therm.xml
\brief structure 3D orthotrope (2 matériaux) avec chargement en température uniquement
*/

/**\example data_Qstat_iso.xml
\brief structure 3D isotrope en quasistatique avec évolution de l'effort imposé
*/

/**\example data_Qstat_iso_contact.xml
\brief structure 3D isotrope avec des interfaces de type contact avec frottement et sans jeu en quasitatique avec évolution du chargement
*/

/**\example data_Qstat_iso_jeu.xml
\brief structure 3D isotrope avec des interfaces de type contact avec frottement et jeu en quasitatique avec évolution du chargement
*/

