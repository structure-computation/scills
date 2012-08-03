#ifdef METIL_COMP_DIRECTIVE
    #pragma src_file multiscale_geometry_mesh.cpp
    #pragma src_file assignation_material_properties_Sst.cpp
    #pragma src_file assignation_material_properties_Interface.cpp
    #pragma src_file multiscale_operateurs.cpp
    #if DIM == 2
        #pragma src_file formulation_2_double_elasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_elasticity_orthotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_elasticity_damageable_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_plasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_2_double_mesomodele.cpp
    #else
        #pragma src_file formulation_3_double_elasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_elasticity_orthotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_elasticity_damageable_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_plasticity_isotropy_stat_Qstat.cpp
        #pragma src_file formulation_3_double_mesomodele.cpp
    #endif
    #pragma src_file iterate_stat_Qstat.cpp
    #pragma src_file affichage.cpp
#endif



//lecture des parametres
#include "DEFINITIONS/Process.h"

using namespace Metil;

#ifndef INFO_TIME
#define INFO_TIME
#endif 

#include "../LMT/include/io/ioexception.h"
using namespace LMT;
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
    /// Verification des arguments
    if ( argc != 3 and argc != 4 ) {
        std::cerr << "usage : nom_executable + dimension +  nom_complet_du_fichier_de_donnees. (+ mpi)" << std::endl;
        std::cerr << "ex : ./SC_multi_2 num_model num_calcul  (+ mpi)" << std::endl;
        return 1;
    }
    try {
        Sc2String id_model = argv[1];
        Sc2String id_calcul = argv[2];
        
        Process process;
        process.initialisation_MPI(argc, argv);
        process.lecture_fichiers(id_model, id_calcul);
        

        if(process.data_user->options.mode=="visu_CL"){
            if(process.parallelisation->is_master_cpu()) std::cout << "Mode de visualisation des bords" << std::endl;
            process.geometry_user->visualize_group_edges_within_geometry(*process.data_user);
        }
        else{
            process.geometry_user->split_group_edges_within_geometry(*process.data_user);
            process.preparation_calcul();
            process.boucle_multi_resolution();
        }
        process.finalisation_MPI();
        if(process.parallelisation->is_master_cpu()) {
            std::cout << "End of SC_multi_" << DIM << ".exe " << id_model << " " << id_calcul  << std::endl;
        }
    } catch ( const IoException &e ) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}




/** \mainpage Fichier general pour le programme multiechelle.
 
Ce programme est constitue de differentes sections chacune correspondant à un fichier cpp separe :
 
- \ref lecture "Lecture des donnees"
 
- \ref Maillage_geometrie "Creation des sous-structures et interfaces"
 
- \ref Materiaux "Assignation des materiaux aux sous-structures et des proprietes aux interfaces"
 
- \ref Operateurs "Creation des operateurs"
 
- \ref Strategie_iterative "Strategie iterative"
 
- \ref Exemples "Exemples de fichiers de donnees xml"
*/


/** \page lecture Lecture des donnees.
 
Dans un premier temps le fichier xml entre comme argument est lu et 
les differentes quantites sont affectees aux instances des classes correspondantes. 
 
L'ordre dans lequel sont rentrees les differents champs xml n'a pas d'importance.
 
Les champs du fichier xml doivent être inclus dans  :
 
\code  
<?xml version="1.0"?>
<root>
 
</root>
\endcode
 
\section Lecture_calcul Lecture des donnees liees au calcul.
La lecture du xml permet d'assigner les differents champs de l'instance process de la class Process.
Les champs xml suivants sont lus (cliquer sur les liens pour voir le detail des sous-champs à remplir):
- paramètres assignes à la classe MultiScaleData et Process \ref parametres_xml    
\code <parametres /> \endcode               
- paramètres temporels assignes à la classe TimeData. \ref parametres_temporels_xml
\code <parametres_temporels /> \endcode  
- parametres d'affichage (classe SavingData). \ref parametres_affichage_xml
\code <parametres_affichage /> \endcode  
- paramètres pour les directions de recherche (classe MultiScaleData). \ref parametres_recherche_xml 
\code <parametres_recherche /> \endcode  
 
\section Lecture_structure Lecture des donnees associees à la geometrie et au maillage 
On lit ensuite la ligne suivante correspondant aux maillages et nombre de maillages à lire que l'on assigne à l'instance de la classe GeneralData contenue dans Process :\ref mesh
\code 
    <mesh />
\endcode
 
\section Lecture_CL Lecture des conditions aux limites.
On termine par la lecture des differentes lignes de conditions aux limites (1 bloc CL par condition ce qui correspond à une entree dans le vecteur de conditions aux limites classe Boundary) : \ref CL
\code 
   <CL>
   
   </CL>
\endcode
 
 
 
 
\section Details Details des champs xml
 
\subsection parametres_xml - Details parametres
Liste exhaustive des champs disponibles :
\code <parametres sous_integration="1" type_sous_integration="p" nbmacro="9" multiechelle="1" blocage_modes_rigides="0" mvts_bloques="Ty Tx Rz" nbitermax="800" facteur_relaxation="0.8" critere_erreur="1e-5" type_erreur="ddr" nb_threads="1" calcul_micro="non" deltaT="0"/>   \endcode
- sous-integration (obligatoire) : 1 pour utiliser la sous-integration 0 sinon.
- type_sous_integration (obligatoire) : dans le cas ou sous-integration=="oui" on lit ce paramètre qui prend la valeur "h" pour un redecoupage des elements (ne marche que pour les triangles quad et tetra) ou la valeur "p" pour une augmentation du degre d'interpolation de l'element (passage de P1 à P2 par ajout de noeuds partout)
- nbmacro (obligatoire) : nombre de fonctions macro spatiales pour toutes les interfaces (9 en 3D et 4 en 2D).
- multiechelle (obligatoire) : booleen = 1 pour effectuer un calcul multiechelle, 0 pour un calcul monoechelle
- blocage_modes_rigides (obligatoire) : booleen =1 pour dire de bloquer explicitement les modes de corps rigides.
- mvts_bloques (obligatoire) : vecteur contenant les valeurs Tx Ty Tz Rx Ry Rz (separes par des espaces). A renseigner uniquement si blocage_modes_rigides == 1. Par defaut lorsqu'il n'y a que des interfaces de type effort impose sur le bord, le programme determine automatiquement les mouvements macro à bloquer. Par contre s'il ne faut bloquer qu'une translation (exemple lorsqu'on utilise une condition de symetrie), il est necessaire d'aider le programme en indiquant la translation à bloquer.
- nbitermax (obligatoire) : nombre d'iterations maximal pour la LATIN. Si le critère est atteint avant on sort de la boucle sur les iterations.
- facteur_relaxation (obligatoire) : 0.8 est souvent une bonne valeur
- critere_erreur (obligatoire) : critère à atteindre pour arreter l'algorithme latin.
- type d'erreur (obligatoire) : choix possible entre "ddr" et "energie" (cf. \ref calcul_erreur ).
- nb_threads (obligatoire) : nombre de processus à lancer sur une même machine
- deltaT (obligatoire) : valeur du champ de temperature pour un calcul thermoelastique
 
\subsection parametres_temporels_xml - Details parametres_temporels
Liste exhaustive des champs disponibles :
\code <parametres_temporels type_de_calcul="Qstat" nbpastemps="24" pasdetemps="0.05" /> \endcode
- type_de_calcul (obligatoire) : choix entre "Qstat" pour la quasistatique et "stat" pour la statique.
- nbpastemps (non necessaire en statique) : nombre d'intervalle d'etude (ou piquet de temps)
- pasdetemps (non necessaire en statique) : valeur du pas de temps
 
\subsection parametres_affichage_xml - Details parametres_affichage
Liste exhaustive des champs disponibles :
\code   <parametres_affichage interactivite="1" affich_mesh="0" affich_resultat="1" type_affichage="Sinterieur" save_or_display="display" display_fields="dep qtrans num sigma epsilon ener" list_error="1" display_error="1" repertoire_save="./tmp/3d_2_cubes/" name_data="data1" affich_depl_pt="1" coor_point="100 100 100"/> \endcode 
- interactivite (obligatoire) : booleen = 1 permettant à la fin d'un calcul de post-traiter les resultats (cf. Interactivite), 0 sinon et par defaut on utilise alors les autres champs pour afficher le resultat ou non.
- affich_mesh (obligatoire) : (booleen =1 ) une fois les maillages crees et surdiscretises pour les ssts et les interfaces, il est possible de visualiser les maillages avec eclate des sous-structures ou interfaces, type de materiau ou d'interface et numero.
- type_affichage (obligatoire) : choix parmi "Sinterieur", "Sbord" et "Inter" pour afficher le maillage des sous-structures ou de bord des ssts ou des interfaces.
- save_or_display (obligatoire) : choix parmi "save" ou "display" : save pour ne sortir que le fichier paraview0.vtu ou display pour le creer et lancer paraview.
- affich_resultat (obligatoire) : (booleen = 1) affichage du resultat après critère d'erreur atteint. On cree alors un fichier pvd permettant de visualiser la deformee ou autre au cours des pas de temps.
- display_fields (obligatoire) : champs separes par des espaces permettant de selectionner les champs à afficher. Pour la liste des noms disponibles, se reporter aux fichiers de formulations dans le repertoire FORMULATIONS. 
- repertoire_save (obligatoire) : repertoire de sauvegarde des resultats de calcul.
- name_data (obligatoire) : nom generique pour chaque pas de calcul (ex : pour le pas 2 on aura repertoire_save/data1_2.vtu). Le fichier pvd cree sera lui aussi dans le repertoire sous le nom : repertoire_save/data1.pvd
- list_error (obligatoire) : booleen = 1 pour faire defiler l'erreur au cours des iterations
- display_error (obligatoire) : booleen = 1 pour afficher la courbe d'erreur en fin de calcul
 
\subsection parametres_recherche_xml - Details parametres_recherche
Liste exhaustive des champs disponibles :
\code  <direction_recherche ktype="scalaire_auto_CL" kfact="1" formule="E/L" copydirection="1" /> \endcode
- ktype (obligatoire) : type de direction de recherche : les noms disponibles sont detailles dans LATIN::ktype, le mieux est d'utiliser "scalaire_auto_CL" avec kfact=1
- kfact (obligatoire) : facteur multiplicatif des directions de recherche ou directement la direction de recherche
- copydirection (obligatoire) : booleen =1 : on effectue une moyenne des directions obtenues que l'on assigne aux deux cotes de l'interface
 
\subsection mesh - Details mesh
Liste exhaustive des champs disponibles :
\code  <mesh repertoire_des_maillages="./DATAFILES/test3D/3d_isotrope_cubes/" nom_fichier_qualification_materiaux="Vqualif.py" nom_des_maillages="V" nb_de_maillages="2" extension=".avs" jeu_physique="0" /> \endcode
- repertoire_des_maillages (obligatoire) : contient les maillages des sous-structures à lire
- nom_des_maillages (obligatoire) : nom generique pour les maillages (ex : V0.avs, V1.avs ...). Le premier numero est necessairement 0.
- extension (obligatoire) : permet de lire des fichiers differents selon l'extension
- nb_de_maillages (obligatoire) : correspond au nombre de sous-structures
- jeu_physique (obligatoire) : booleen = 1 pour lire les maillages d'interfaces disjointes
- nom_maillages_jeu (uniquement si jeu_physique==1) : nom generique des maillages d'interfaces donnes. Les noms doivent être generes de la facon suivantes : nom_maillages_jeu0-0-2.avs -> interface extraite de la sst numero 0 et ayant comme voisins les sst 0 et 2.
- inter_jeu (uniquement si jeu_physique==1) : vecteur de paire de sous-structures separes par des ; (ex : "0 2;1 3")
 
 \subsection CL - Details CL
 
 A l'interieur d'un bloc plusieurs champs sont disponibles selon le type de CL
\code     
<CL>
      <parametres box="100 0 100    100 100 200" unit="mm" comp="effort" />
      <fct_spatiale fonction="1;0;0" />
      <fct_temporelle intervalle="0. 0.1"  fonction="0" />
      <fct_temporelle intervalle="0.1 1.1"  fonction="-50*sin(3.14159265358979*(t-0.1))" />
   </CL>
 \endcode 
   - le champ box contient une boite definie par ses deux extremites (2*2 composantes en 2d et 3*2 en 3d) permettant d'extraire le maillage inclue dans cette boite
   - comp renseigne sur le type de comportement (effort, depl, sym, depl_normal)
   - le champ fct_spatiale n'est utilise que pour les interfaces de type effort ou deplacement. A savoir, on rentre les valeurs des fonctions imposees selon x, y (et z en 3d) separees par des ; (ne pas mettre d'espace). Les fonctions doivent être des fonctions de x, y (et z en 3d) selon la zone consideree. Ex : surface (0 0  0 1 ) ->  fonction de y
   Pour une condition de type deplacement normal impose (depl_normal) seule une fonction est necessaire, on doit aussi ajouter le champ normale (x,y ou z). Pour une condition de symetrie ce champ peut être efface.
   \code 
   <fct_spatiale fonction="(0*y/2+0)" normale="x"/>
   \endcode
   - On renseigne autant de blocs fct_temporelle que d'intervalles pour lesquels la fonction temporelle est differente. La valeur de la fonction est evaluee pour chaque pas de temps et multipliee à la fonction spatiale.
   
   Conseil : Il est preferable de prendre une fct_spatiale unitaire pour que le produit fct_spatiale*fct_temporelle soit simple à interpreter.
 
 */



/** \page Exemples Exemples de fichiers de donnees xml
 
Voici une liste d'exemples de fichiers de donnees xml compris dans le repertoire DATAFILES
 
Le premier exemple est une structure 2D chargee uniquement par des efforts ou deplacements (la variation de temperature est nulle). On presente ici deux possibilites : 
 - contrainte plane avec chargement en effort :  \ref data_stat_isotrope_CP_effort.xml
 - deformation plane avec chargement normal (deformation longitudinale) : \ref data_stat_isotrope_DP_depl.xml
 
Le deuxième exemple est une structure 3D orthotrope (2 materiaux) avec chargement en temperature uniquement : \ref data_stat_ortho_therm.xml
 
Le troisième exemple est une structure 3D isotrope en quasistatique avec evolution de l'effort impose : \ref data_Qstat_iso.xml
 
Le quatrième exemple est une structure 3D isotrope avec des interfaces de type contact avec frottement et sans jeu en quasitatique avec evolution du chargement : \ref data_Qstat_iso_contact.xml
 
Le dernier exemple est une structure 3D isotrope avec des interfaces de type jeu (contact avec frottement) en quasistatique avec evolution du chargement :  \ref data_Qstat_iso_jeu.xml
 
*/

/** \example data_stat_isotrope_CP_effort.xml
\brief Structure 2D chargee uniquement par des efforts ou deplacements (la variation de temperature est nulle) : contrainte plane avec chargement en effort
*/

/** \example data_stat_isotrope_DP_depl.xml
\brief Structure 2D chargee uniquement par des deplacements normaux (la variation de temperature est nulle) : deformation plane
*/

/** \example data_stat_ortho_therm.xml
\brief structure 3D orthotrope (2 materiaux) avec chargement en temperature uniquement
*/

/**\example data_Qstat_iso.xml
\brief structure 3D isotrope en quasistatique avec evolution de l'effort impose
*/

/**\example data_Qstat_iso_contact.xml
\brief structure 3D isotrope avec des interfaces de type contact avec frottement et sans jeu en quasitatique avec evolution du chargement
*/

/**\example data_Qstat_iso_jeu.xml
\brief structure 3D isotrope avec des interfaces de type contact avec frottement et jeu en quasitatique avec evolution du chargement
*/

