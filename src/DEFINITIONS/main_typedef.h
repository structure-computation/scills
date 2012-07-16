#ifndef LMT_INCLUDES_H
#define LMT_INCLUDES_H

/// Definition des types basiques
#include "../UTILS/Sc2String.h"
typedef int Id;
//typedef Sc2String String;   /// Le type des chaines de caracteres Incompatible avec Metil
typedef TYPEREEL Scalar;

/// Definition des types de vecteurs
#include "../../LMT/include/containers/vec.h"
#include "../../LMT/include/containers/vecpointedvalues.h"
typedef LMT::Vec<Scalar,DIM> Point;                     /// Defini un point (ou un vecteur) de l'espace R2 ou R3 (selon la dimension du probleme)
typedef LMT::Vec<Scalar> Vector;                        /// Un vecteur de reels de taille dynamique
typedef LMT::Vec<Scalar,DIM*(DIM+1)/2> VoigtVector;     /// Un vecteur dans les notations de Voigt

/// Definition des type de matrices
#include "../../LMT/include/containers/mat.h"
typedef LMT::Mat<Scalar,LMT::Gen<>,LMT::Dense<> > DenseMatrix;                      /// Une matrice carree dense de reels et de taille dynamique
typedef LMT::Mat<Scalar,LMT::Gen<>,LMT::SparseLine<> > SparseMatrix;                /// Une matrice carree bande de reels et de taille dynamique
typedef LMT::Mat<Scalar,LMT::Sym<>,LMT::SparseLine<> > SymetricMatrix;              /// Une matrice bande symetrique de reels et de taille dynamique
typedef LMT::Mat<Scalar,LMT::Sym<DIM*(DIM+1)/2>,LMT::SparseLine<> > VoigtMatrix;    /// Une matrice dans les notations de Voigt
typedef LMT::Mat<bool,LMT::Gen<>,LMT::SparseLine<> > SymbolicMatrix;                /// Une matrice carree symbolique (de booleens) de taille dynamique

/// Definition des types de maillage
#include "meshmulti.h"
#include "../MAILLAGE/meshcaracinter.h"
#include "../../build/problem_pb_elast/mesh_carac.h"
#include "../../build/problem_pb_elast/problem.h"
#include "../../LMT/include/formulation/formulation_ancestor.h"
typedef LMT::Problem_pb_elast<Scalar,DIM> Problem;
typedef LMT::FormulationAncestor<Scalar> Formulation;
typedef Meshmulti<LMT::Mesh_carac_pb_elast<Scalar,DIM> > SstMesh;   /// type de maillage pour les sous-structures
typedef LMT::Mesh<LMT::Meshcaracinter<DIM,DIM-1> > EdgeMesh;        /// type de maillage pour les bords des sous-structures
typedef LMT::Mesh<LMT::Meshcaracinter<DIM,DIM-1>,0 > InterfaceMesh; /// Type de maillage pour les interfaces


#endif  //LMT_INCLUDES_H