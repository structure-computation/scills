// fonctions speciales
#include "algebre.h"
#include "utilitaires.h"

// fonction utilisees pour le maillage et la geometrie
#include "calculate_measure_G_SST.h"
#include "create_mesure_G_neq_INTER.h"
#include "create_SST_INTER.h"
#include "p_surdiscretisation.h"
#include "surdiscretisation.h"
#include "p_create_SST_edge.h"

// #include "affichage.h"
using namespace LMT;
using namespace std;

//***************************************************
// Programme principal : creation maillage geometrie
//***************************************************

/** \ingroup Maillage_geometrie
\brief Fonction principale de création des maillages et géométries
 
- Dans un premier temps on crée les sous-structures et interfaces à partir des informations contenues dans la class STRUCTURE (inclus dans Param) et des conditions aux limites CL. On obtient ainsi un vecteur de sous-structures \ref Sous_structures et d'interfaces \ref Interfaces.
La creation des sous-structures et interfaces est effectuée de la facon indiquée dans la fonction create_SST_INTER() .
 
On crée ensuite les données géométriques associées aux interfaces : calculate_measure_G_neq_INTER.
 
- Dans un second temps, on modifie le maillage des sous-structures pour prendre en compte ou non (d'après le paramètre Param::sousint) la sous-integration (cf. surdiscretise_SST ).
 
- Enfin on détermine certaines caractéristiques géométriques des Sst (Sst::G, Sst::measure, Sst::Edge::G).
 
Rq : La réorganisation des noeuds des Sous-structures n'est pas nécessaire étant donné qu'on utilise une renumérotation dans le solveur sparse (symamd).
 
*/
template<class TV1, class TV2, class TV3, class TV4,class TV5> void multiscale_geometry_mesh(const XmlNode &n,TV1 &S, TV2 &Inter, Param &process, TV5 &CL, TV3 &Stot,TV3 &SubS, TV4 &SubI) {

   // creation SST et INTERFACE - maillage
   create_SST_INTER(*process.structure,S,Inter,CL,process,Stot,SubS,SubI);
   
   //creation des normales, centre de gravite et mesure des interfaces
   if (process.rank == 0) cout << "\t Centre de gravite des interfaces + mesure + normales equivalentes" << endl;
   apply_mt(SubI,process.nb_threads,calculate_measure_G_neq_INTER(),S);

   if (process.sousint==1) {
      if (process.rank == 0) cout << "Modification des SST " << endl;
      if (process.rank == 0) cout << "\t Sur-discretisation des elements" << endl;
      if (process.type_sousint=="h") {
         // surdiscretisation des ssts
         apply_mt(SubS,process.nb_threads,surdiscretise_SST());
         // correspondance interface - cote SST - creation des nodeeq
         apply_mt(SubI,process.nb_threads,decoup_INTER(),S);
      } else if(process.type_sousint=="p") {
         //apply_mt(S,process.nb_threads,p_surdiscretise_SST());
          for( unsigned i=0;i<SubS.size() ;i++ ){
//               S[i].mesh.sousint(Number< TV1::template SubType<0>::T::dim-1 >() );
              SubS[i].mesh.sousintegration = "p";
          }
          apply_mt(SubI,process.nb_threads,p_decoup_INTER(),S);
      } else {if (process.rank == 0) cout << "Mauvais choix de sous integration : h ou p " << endl; assert(0);}
      
   } else {
      // correspondance interface - cote SST - creation des nodeeq
      apply_mt(SubI,process.nb_threads,copy_INTER(),S);
   }

   if (process.rank == 0) cout << "\t Centre de gravite des SST et cotes + mesure SST" << endl;
   apply_mt(SubS,process.nb_threads,calculate_measure_G_SST());
   
   //mettre à jour la taille des sst sur tous les pros
   if (process.size > 1 and process.sousint==1)
       for (unsigned i=0;i<S.size();i++){
       for( unsigned k1=0;k1 <process.multi_mpi->repartition_sst.size()  ;k1++ )
           if (find(process.multi_mpi->repartition_sst[k1],_1==(int)i)) {
           MPI_Bcast(&S[i].mesh.node_list_size,1,MPI_INT,k1,MPI_COMM_WORLD);
           MPI_Bcast(&S[i].mesh.elem_list_size,1,MPI_INT,k1,MPI_COMM_WORLD);
           MPI_Bcast(S[i].G.ptr(),S[i].G.size(),MPI_DOUBLE,k1,MPI_COMM_WORLD);
           break;
           }
       }

   if (process.rank == 0)
       calcul_taille_probleme(S, Inter);

}

