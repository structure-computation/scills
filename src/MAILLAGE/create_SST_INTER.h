#ifndef CREATE_SST_INTER
#define CREATE_SST_INTER

#include "../DEFINITIONS/Sst.h"
#include "../DEFINITIONS/Interface.h"
#include "../DEFINITIONS/MacroProblem.h"
#include "../DEFINITIONS/Boundary.h"
#include "Process.h"

#include "../GEOMETRY/GeometryUser.h"
#include "../COMPUTE/DataUser.h"

using namespace Metil;
#include "mpi.h"

extern "C" {
    // #include "metis.h"
//     void METIS_PartGraphRecursive(int *, long long int *, long long int *, long long int *, long long int *, int *, int *, int *, int *, int *, long long int *);
    void METIS_PartGraphRecursive(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
}

using namespace LMT;

/** \defgroup Maillage_geometrie_sst G�om�trie et maillage des Sous-structures
\ingroup Maillage_geometrie
*/
/** \defgroup Maillage_geometrie_inter G�om�trie et maillage des Interfaces
\ingroup Maillage_geometrie
*/


/** \ingroup Maillage_geometrie_sst
\brief Cr�ation du nombre de Sst et affectation du numero du materiau par l'interm�diaire du data_user (fichier json)
 
La structure data_user contient une liste de group_elements correspondant � chaque sst et une id_material. Cet id_material correspond � l'identificateur des mat�riaux donn�s dans le fichier json.
A cette etape, on sp�cifie la taille du vecteur de Sous-structures par l'interm�diaire geometry_user.nb_group_elements et on assigne le num�ro du mat�riau pour chaque sous-structure.
*/
void create_SST_typmat(DataUser &data_user, GeometryUser &geometry_user,Substructures &S,Process &process);


/**\ingroup Maillage_geometrie_sst
\brief Lecture et assignation d'un maillage par sst
 
Le maillage de chaque sst est lu � partir des donn�es charg�e dans geometry_user. Ces donn�es viennent de SC_create_2
 
*/
///Fonction generique
template<class TE, class TM, class TR>
void add_new_elem(TE &e, TM &m, TR &rep_nodes) {}
///Fonctions specialisees
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Tetra,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Tetra(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Wedge,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Wedge(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Hexa,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Hexa(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Triangle,TNB,TN,TD,NET> &e, TM &m, TR&rep_nodes) {
    m.add_element(Triangle(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Quad,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Quad(),DefaultBehavior(),rep_nodes.ptr() );
}
template<class TNB,class TN,class TD,unsigned NET, class TM, class TR>
void add_new_elem(Element<Bar,TNB,TN,TD,NET> &e, TM &m, TR &rep_nodes) {
    m.add_element(Bar(),DefaultBehavior(),rep_nodes.ptr() );
}

//Structure utilisee pour realiser une operation sur les elements du maillage de peau d'une sous-structure : 
//      copie des connectivites des elements de peau du maillage local d'une sous-structure (LMT::mesh) vers l'objet geometry_user (GPU)
struct ConvertMeshConnectivitiesSkin{
   template<class TE> void operator()(TE &e, GeometryUser &geometry_user, int id_sst) const{
      for(unsigned i=0;i<e.nb_nodes;i++){
          geometry_user.find_group_elements(id_sst)->local_connectivities_skin[i][e.number]=e.node(i)->number_in_original_mesh();
      }
   }
};

//conversion du maillage de peau d'une sous-structure (LMT::mesh) aux group_elements
void convert_mesh_skin_to_geometry_user(Sst &S, GeometryUser &geometry_user);


void create_maillage_SST(DataUser &data_user, GeometryUser &geometry_user, Substructures &S, Process &process);


/** \ingroup Maillage_geometrie_inter
\brief Cr�ation des interfaces comprises entre les Sst. Par d�faut elles sont consid�r�es comme parfaites.
 
Le maillage de l'interface est obtenu par la lecture de geometry_user issue de SC_create_2
 
*/

void read_mesh_interface_geometry_user(SstMesh::TM &mesh, GeometryUser &geometry_user, int num_inter) throw(std::runtime_error);

void create_perfect_interfaces(DataUser &data_user, GeometryUser &geometry_user, Substructures &S, VecInterfaces &Inter, Process &process);

/** \ingroup Maillage_geometrie_inter
\brief Cr�ation des interfaces contenues dans une condition aux limites donn�e
 
On recherche les Ssts ayant une intersection avec chacune des conditions aux limites (intersection de boites). On extrait alors du maillage de peau de chaque Ssts s�lectionn�e le maillage contenu dans une condition aux limites donn�e. 
*/
void create_interfaces_CL(DataUser &data_user, GeometryUser &geometry_user, Substructures &S, VecInterfaces &Inter, Boundaries &CL, Process &process);

struct make_interface_inter {
    void operator()(Interface &SubI, GeometryUser &geometry_user) const;
};

struct make_interface_CL {
    void operator()(Interface &SubI, GeometryUser &geometry_user) const;
};

struct mesh_unload {
    void operator()(Sst &S) const;
};

//******************************************************************************************************
// Programme permettant la construction des SST et Interfaces aveclecture de maillage issu de castem
//******************************************************************************************************
/** \ingroup Maillage_geometrie
\brief Cr�ation des Sous-structures et des Interfaces (maillages).
 
On construit ici le vecteur des sous-structures et des interfaces � partir des informations du fichier xml.
 
- Assignation du type de mat�riau aux sous-structures create_SST_typmat().
- Lecture du maillage des sous-structures create_maillage_SST().
- creation des interfaces parfaites create_perfect_interfaces().
- creation des interfaces avec jeu physique create_gap_interfaces().
- creation des interfaces de conditions aux limites. create_interfaces_CL().
*/
                     
#ifndef INFO_TIME
#define INFO_TIME
#endif 
#include "containers/evaluate_nb_cycles.h"

void create_SST_INTER(DataUser              &data_user,
                      GeometryUser          &geometry_user,
                      Substructures         &S,
                      VecInterfaces         &Inter,
                      Boundaries            &CL,
                      Process               &process,
                      PointedSubstructures  &Stot,
                      PointedSubstructures  &SubS,
                      PointedInterfaces     &SubI);


#endif //CREATE_SST_INTER