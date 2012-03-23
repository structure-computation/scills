#ifndef DEFINITION_DECLARATIONS_H
#define DEFINITION_DECLARATIONS_H

// SC_multi_2
struct Param;     // container pour les structures de donnees suivantes:
struct AFFICHAGE; // parametres d'affichage
struct STRUCTURE; // parametres de maillages et geometries
struct LATIN;     // parametres associes a la resolution LATIN
struct MULTI;     // parametres multiechelles
struct TEMPS;     // parametres temporels
struct PROPERTY;  // parametres de proprietes materielles globales
struct MULTI_MPI; // parametres de MPI

// Structures de donnees d'origine
struct Glob;
struct Sst;
struct SstCarac;
struct Interface;
struct InterCarac;
struct Boundary;
struct PARAM_DAMAGE_INTER;
struct PARAM_COMP_INTER;


#endif //DEFINITION_DECLARATIONS_H