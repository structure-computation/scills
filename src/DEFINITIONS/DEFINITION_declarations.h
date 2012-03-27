#ifndef DEFINITION_DECLARATIONS_H
#define DEFINITION_DECLARATIONS_H

// SC_multi_2
struct Process;     // container pour les structures de donnees suivantes:
struct SaveParameters; // parametres d'affichage
struct GeneralParameters; // parametres de maillages et geometries
struct LatinParameters;     // parametres associes a la resolution LATIN
struct MultiScaleParameters;     // parametres multiechelles
struct TimeParameters;     // parametres temporels
struct PROPERTY;  // parametres de proprietes materielles globales
struct MPIParameters; // parametres de MPI

// Structures de donnees d'origine
struct MacroProblem;
struct Sst;
struct SstCarac;
struct Interface;
struct InterCarac;
struct Boundary;
struct PARAM_DAMAGE_INTER;
struct PARAM_COMP_INTER;


#endif //DEFINITION_DECLARATIONS_H