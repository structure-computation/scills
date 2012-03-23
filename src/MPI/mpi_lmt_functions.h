//
// C++ Interface: mpi_lmt_functions
//
// Description: Fonction MPI pour passer les classes LMT 
//
//
// Author: Alain Caignot <caignot@lmt.ens-cachan.fr>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MPI_LMT_H
#define MPI_LMT_H

#include "mpi.h"

/** \example send_recv.cpp Exemple de passage d'�l�ments de la classe LMT d'un ordinateur � un autre...

Dans cet exemple, nous pr�sentons comment on utilise les fonctions MPI_Send et MPI_Recv sur des donn�es �l�mentaires, sur les vecteurs, sur les vecteurs de vecteurs, sur des matrices full et sparse, et enfin sur la classe Time du code COFAST+ de D. Violeau.
 
 
*/
/** \example somme_partiel_series.cpp Exemple du calcul d'une somme partielle de s�rie sur plusieurs processeurs.

Dans cet exemple, on traite le calcul de la somme partiel d'une suite sur plusieurs processeurs. Cet exemple permet de traiter la synchronisation par 2 moyens diff�rents, il pr�sente de nouvelles fonctions MPI standard.
*/


/** \defgroup MPI Parall�lisation avec MPI et la biblioth�que LMT
\brief Sp�cifications des commandes MPI usuelles appliqu�es � la biblioth�que LMT

La syntaxe des fonctions MPI �tant un peu lourde, pour simplifier son utilisation dans les cas usuels que nous pouvons rencontrer, nous avons surd�finie ces fonctions. Nous vous pr�sentons les diff�rentes fonctions et leurs utilisations class�s par type : envoie, r�ception, ...
Seules les fonctions MPI de base ont �t� cod�e � savoir l'envoie avec \ref MPI_Send et  la r�ception avec \ref MPI_Recv. Mais il est tout � fait possible de recr��er de nouvelles fonctions en suivant le formalisme adopt� et en se basant sur les fonctions d�j� d�finies.

Les objets g�r�s pour le moment sont les donn�es �l�mentaires (\c int, \c double, ...), les vecteurs \c LMT::Vec, les matrices \c LMT::Mat \c full ou \c sparse. Il y a �galement un exemple avec la structure \c Time du code Multi de D.Violeau. 

Ne sont donc pas g�r�s les matrices � stockage \c skyline (qui ne sont pas utilis�s dans le code Multi). De mani�re g�n�ral, tous les conteneurs peuvent �tre envoy�s plus ou moins facilement. Ainsi que des structures homog�nes (contenant que des �l�ments du m�me type � l'int�rieure ou de type pouvant �tre reconstruit, ie. passer tout en \c double et transformer certaines parties du r�sultats en \c int). En revanche, il n'est pas raisonnable de vouloir passer des structures fortement h�t�rog�nes comme la classe \c LMT::Mesh !
*/

/**
\mainpage Documentation des fonctions MPI adapt�s � la biblioth�que LMT

\section intro Introduction
La n�cessit� de parall�liser les codes de calcul se fait ressentir lorsque l'on veut passer des gros calculs et ceci dans tous les domaines. Il y a deux fa�ons de parrall�liser un code en fonction de l'architecture mat�riel que l'on souhaite utiliser : un machine parall�le � m�moire partager ou un cluster. Pour ce qui est des machines SMP, la classe LMT permet facilement de parall�liser son code en rempla�ant les boucles \c for par la fonction \c apply_mt.

Le standard MPI (Message Passing Interface) permet de parall�liser le code sur des clusters de machines, ensemble qui sont de 5 � 10 fois moins co�teux que les machines SMP.

Nous avons donc mis en place des fonctions de base permettant de simplifier l'appel aux fonctions MPI standard avec les �l�ments standard de la biblioth�que LMT.

\section Rappel Rappel des commandes �l�mentaires � connaitre pour pouvoir utiliser LAM/MPI :

\subsection Installation Installation : 
- il faut installer les fichiers de d�veloppement : librairies, sources, compilateur (<tt>mpi.h, mpic++, mpicc, ...</tt>)
- il faut installer les utilitaires qui permettent de lancer les processus MPI (<tt>lamboot, lamwipe, lamclean, ...</tt>)

\subsection Executer Lancer un programme :
- dans le fichier \c lamhosts, mettez la liste des machines, la machine locale devant �tre mise en premier, elle sera consid�r� comme le master
- <tt>lamboot -v lamhosts</tt> permet d'executer le serveur MPI sur toutes les machines
- <tt>mpirun -v -np nbr_proc ./programme arguments</tt> permet de lancer le programme
- <tt>lamclean -v lamhosts</tt>, et <tt>lamwipe -v lamhosts</tt> permet de fermer toutes les instances MPI � la fin du programme pour nettoyer les machines

\subsection Developper D�velopper : 
- un programme doit commencer par :
\code
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size(MPI_COMM_WORLD, &size);
\endcode
ce qui permet d'initialiser MPI et � chaque noeuds de connaitre son environnement (nombre de processus, son num�ro)

- tous les programmes doivent se terminer par : 
\code
    MPI_Finalize();
    return 0;
\endcode

- les commandes permettant d'envoyer des messages : \c MPI_Send, \c MPI_Isend
- les commandes permettant de recevoir des messages : \c MPI_Send, \c MPI_Isend
- les commandes permettant de savoir si un message est l� : \c MPI_Probe, \c MPI_IProbe
- les commandes de synchronisation : \c MPI_Barrier
- ...

Plus de d�tails sur MPI http://www.lam-mpi.org/

\section MPI_LMT Liens vers les sections documentant les fonctions MPI pour la biblioth�que LMT
\ref MPI 
 - D�finition des types : \ref Mpitype
 - Envoyer des donn�es : \ref MPI_Send
 - Recevoir des donn�es : \ref MPI_Recv
 - Apprendre par l'exemple : <a class="el" href="examples.html">Exemples</a>

*/

using namespace LMT;



/** \ingroup MPI
\defgroup Mpitype Mpi_type
\brief Description des types de base MPI.

MPI utilise des types standard. Cette classe permet de traduire un type �l�mentaire en type MPI. 
Les types MPI standard disponibles sont double, int, char, short, long, float, long double, unsigned, unsigned char, unsigned short, unsigned long.
Les types MPI_PACKED, MPI_BITES, MPI_BOOL ne sont pas support�s.
*/
/**
\ingroup Mpitype

MPI utilise des types standard. Cette classe permet de traduire un type �l�mentaire en type MPI. 
Les types MPI standard disponibles sont double, int, char, short, long, float, long double, unsigned, unsigned char, unsigned short, unsigned long.
*/
template<class T> struct Mpi_type {};
template<> struct Mpi_type<double>            { static MPI_Datatype res() { return  MPI_DOUBLE        ;}};
template<> struct Mpi_type<int>               { static MPI_Datatype res() { return  MPI_INT           ;}};
template<> struct Mpi_type<char>              { static MPI_Datatype res() { return  MPI_CHAR          ;}};
template<> struct Mpi_type<short int>         { static MPI_Datatype res() { return  MPI_SHORT         ;}};
template<> struct Mpi_type<long int>          { static MPI_Datatype res() { return  MPI_LONG          ;}};
template<> struct Mpi_type<float>             { static MPI_Datatype res() { return  MPI_FLOAT         ;}};
template<> struct Mpi_type<long double>       { static MPI_Datatype res() { return  MPI_LONG_DOUBLE   ;}};
template<> struct Mpi_type<unsigned int>      { static MPI_Datatype res() { return  MPI_UNSIGNED      ;}};
template<> struct Mpi_type<unsigned char>     { static MPI_Datatype res() { return  MPI_UNSIGNED_CHAR ;}};
template<> struct Mpi_type<unsigned short int>{ static MPI_Datatype res() { return  MPI_UNSIGNED_SHORT;}};
template<> struct Mpi_type<unsigned long int> { static MPI_Datatype res() { return  MPI_UNSIGNED_LONG ;}};




/** \ingroup MPI
\defgroup MPI_Send 
\brief Fonctions pour envoy�s des donn�es

Seuls les fonctions MPI_Send de MPI ont �t� surd�finie. Il existe d'autres types de fonctions d'envoie comme les envoies non bloquant MPI_ISend qui n'ont pas �t� cod�e. Il est cependant ais� de les red�finir pour ceux qui en ont besoin.
*/

/** \ingroup MPI_Send
\brief Passer un �l�ment de type standard - support du tag

Il suffit de lancer la commande MPI_Send(�l�ment, rank) pour envoyer l'�l�ment sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type
*/
template<class TN>
void MPI_Send(TN &a,int rank,int tag=201) {
    MPI_Send(&a, 1, Mpi_type<TN>::res() , rank, tag, MPI_COMM_WORLD);
}

/** \ingroup MPI_Send
\brief Passer un vecteur de la classe LMT de type standard - support du tag

Il suffit de lancer la commande MPI_Send(vecteur, rank) pour envoyer le vecteur sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type
*/
template<class TN>
void MPI_Send(Vec<TN> &vec,int rank,int tag=201) {

    MPI_Send(vec.ptr(), vec.size(),
Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD);
}


/** \ingroup MPI_Send
\brief Passer un vecteur de la classe LMT de type standard
 
Il suffit de lancer la commande MPI_Isend(vecteur, rank) pour envoyer le vecteur sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type
*/
template<class TN>
void MPI_Isend(Vec<TN> &vec,int rank,MPI_Request &request, int tag=201) {
    MPI_Isend(vec.ptr(), vec.size(),
              Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD,&request);
}





/** \ingroup MPI_Send
\brief Passer un vecteur de la classe LMT de vecteur de la classe LMT de type standard

Il suffit de lancer la commande MPI_Send(vecteur, rank) pour envoyer le vecteur sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type. Le sous-type des vecteurs doit �tre le m�me, cependant la taille des sous-vecteurs peut �tre diff�rente.
*/
template<class TN>
void MPI_Send(Vec<Vec<TN> > &vec,int rank,int tag=201) {
    Vec<TN> a;
    int size=0;
    for( int i=0;i<vec.size() ;i++ ){
        size+=vec[i].size();
    }
    a.resize(size);size=0;
    for( int i=0;i<vec.size() ;i++ ){
        a[(range(vec[i].size())+ size)]= vec[i];
        size+=vec[i].size();
    }
    MPI_Send(a.ptr(), a.size(),
Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD);
}
/** \ingroup MPI_Send
\brief Passer une matrice de la classe LMT de type standard

Il suffit de lancer la commande MPI_Send(matrice, rank) pour envoyer la matrice sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type.
*/
template<class TN>
void MPI_Send(Mat<TN> &mat,int rank,int tag=201) {
    MPI_Status status;
    MPI_Send(mat.data.ptr(), mat.nb_cols()*mat.nb_rows(),
Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD);
}


/** \ingroup MPI_Send
\brief Passer une matrice de la classe LMT de type sparse - support du tag

Il suffit de lancer la commande MPI_Send(matrice, rank) pour envoyer la matrice sur la machine rank. Tous les types standard qui peuvent �tre pass� sont donn�s dans Mpi_type.
*/
template<class TN, class TO>
void MPI_Send(Mat<TN,TO,SparseLine<Col> > &mat,int rank,int tag=201) {
    MPI_Status status;
    Vec<int> row;row.resize(mat.data.size());
    int size=0;
    for( int i=0 ;i<mat.data.size() ;i++ ){
         row[i]=mat.data[i].indices.size();
         size+=mat.data[i].indices.size();
    }
    Vec<TN> a; a.resize(row.size() + 2*size+2);
    a[0]=mat.nb_rows();a[1]=mat.nb_cols();int pos=2;
    for( int i=0 ;i<mat.data.size() ;i++ ){
         a[pos]=row[i];pos+=1;
         a[range(row[i])+pos]=mat.data[i].indices;pos+=row[i];
         a[range(row[i])+pos]=mat.data[i].data;pos+=row[i];
    }
    MPI_Send(a.ptr(), a.size(),
Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD);
}



/*/** \ingroup MPI_Send
\brief Passer un vecteur de la classe Time

Il suffit de lancer la commande MPI_Send(vecteur, rank, topass) pour envoyer le vecteur sur la machine rank o� topass est un vecteur de string contenant le nom des composantes � faire passer de la class Time.
*/
// template<class TN>
// void MPI_Send(Vec< Time<TN> > &t,int rank, Vec<string> topass){
// MPI_Send(t,rank,topass,201);
// }
/*/** \ingroup MPI_Send
\brief Passer un vecteur de la classe Time - support du tag
*/
// template<class TN>
// void MPI_Send(Vec< Time<TN> > &t,int rank, Vec<string> topass,int tag){
//     Vec<TN> a;
//     a.resize( t[0].W1.size() * topass.size() * t.size());
//     int size=0;
//     for( int j =0 ;j < t.size() ; j++){
//         for( int i=0; i < topass.size() ;i++ ){
//             if (topass[i] == "Wd") {
//                 a[range(t[0].W1.size())+size]=t[j].Wd;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Qd") {
//                 a[range(t[0].W1.size())+size]=t[j].Qd;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wtemp") {
//                 a[range(t[0].W1.size())+size]=t[j].Wtemp;size+=t[0].W1.size();
//             }
//             if (topass[i] == "F") {
//                 a[range(t[0].W1.size())+size]=t[j].F;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wp") {
//                 a[range(t[0].W1.size())+size]=t[j].Wp;size+=t[0].W1.size();
//             }
//             if (topass[i] == "W") {
//                 a[range(t[0].W1.size())+size]=t[j].W;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Fchap") {
//                 a[range(t[0].W1.size())+size]=t[j].Fchap;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wpchap") {
//                 a[range(t[0].W1.size())+size]=t[j].Wpchap;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wchap") {
//                 a[range(t[0].W1.size())+size]=t[j].Wchap;size+=t[0].W1.size();
//             }
//             if (topass[i] == "WtildeM") {
//                 a[range(t[0].W1.size())+size]=t[j].WtildeM;size+=t[0].W1.size();
//             }
//             if (topass[i] == "oldF") {
//                 a[range(t[0].W1.size())+size]=t[j].oldF;size+=t[0].W1.size();
//             }
//             if (topass[i] == "oldWp") {
//                 a[range(t[0].W1.size())+size]=t[j].oldWp;size+=t[0].W1.size();
//             }
//             if (topass[i] == "oldW") {
//                 a[range(t[0].W1.size())+size]=t[j].oldW;size+=t[0].W1.size();
//             }
//             if (topass[i] == "W1") {
//                 a[range(t[0].W1.size())+size]=t[j].W1;size+=t[0].W1.size();
//             }
//             if (topass[i] == "W2") {
//                 a[range(t[0].W1.size())+size]=t[j].W2;size+=t[0].W1.size();
//             }
//             if (topass[i] == "F1") {
//                 a[range(t[0].W1.size())+size]=t[j].F1;size+=t[0].W1.size();
//             }
//             if (topass[i] == "F2") {
//                 a[range(t[0].W1.size())+size]=t[j].F2;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wp1") {
//                 a[range(t[0].W1.size())+size]=t[j].Wp1;size+=t[0].W1.size();
//             }
//             if (topass[i] == "Wp2") {
//                 a[range(t[0].W1.size())+size]=t[j].Wp2;size+=t[0].W1.size();
//             }
//         }
//     }
//     MPI_Send(a.ptr(), a.size(),
//         Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res() , rank, tag, MPI_COMM_WORLD);
// }
















/** \ingroup MPI
\defgroup MPI_Recv
\brief Fonctions pour recevoir des donn�es

Seuls les fonctions MPI_Recv de MPI ont �t� surd�finie. Il existe d'autres types de fonctions d'envoie comme les r�ceptions non bloquant MPI_IRecv qui n'ont pas �t� cod�e. Il est cependant ais� de les red�finir pour ceux qui en ont besoin.
*/
/** \ingroup MPI_Recv
\brief Recevoir un �l�ment standard - support du tag

Il suffit de lancer la commande MPI_Recv(�l�ment, rank) pour recevoir l'�l�ment envoy� dans �l�ment par la machine rank. Il est n�cessaire de connaitre le type � recevoir.
*/
template<class TN>
void MPI_Recv(TN &a,int rank,int tag=201) {
    MPI_Status status;
    MPI_Recv( &a, 1, Mpi_type<TN>::res(), rank, tag, MPI_COMM_WORLD, &status );
}


/** \ingroup MPI_Recv
\brief Recevoir un vecteur de la classe LMT - support du tag

Il suffit de lancer la commande MPI_Recv(vecteur, rank) pour recevoir le vecteur envoy� dans vecteur par la machine rank. Il est n�cessaire de connaitre le type � recevoir. Le vecteur sera resiz� automatiquement.
*/
template<class TN>
void MPI_Recv(Vec<TN> &vec,int rank,int tag=201) {
    MPI_Status status;

    MPI_Probe(rank,MPI_ANY_TAG,MPI_COMM_WORLD, &status );
//    vec.resize( status.st_length / sizeof(TN));//lam-mpi
    vec.resize( status._count / sizeof(TN));//open-mpi
    MPI_Recv( vec.ptr(), vec.size(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &status );
}

/** \ingroup MPI_Recv
\brief Recevoir un vecteur de la classe LMT
 
Il suffit de lancer la commande MPI_Recv(vecteur, rank) pour recevoir le vecteur envoy� dans vecteur par la machine rank. Il est n�cessaire de connaitre le type � recevoir. Le vecteur sera resiz� automatiquement.
*/
template<class TN>
void MPI_Irecv(Vec<TN> &vec,int rank,MPI_Request &request,int tag=201) {
    MPI_Status status;
    MPI_Probe(rank,MPI_ANY_TAG,MPI_COMM_WORLD, &status );
//    vec.resize( status.st_length / sizeof(TN));
    vec.resize( status._count / sizeof(TN));//open-mpi
    MPI_Irecv( vec.ptr(), vec.size(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &request );
}


/** \ingroup MPI_Recv
\brief Recevoir un vecteur de vecteur de la classe LMT - support du tag

Il suffit de lancer la commande MPI_Recv(vecteur, rank) pour recevoir le vecteur de vecteur envoy� dans vecteur par la machine rank. Il est n�cessaire de connaitre le type � recevoir ainsi que les diff�rentes tailles pour le vecteur de reception.
*/
template<class TN>
void MPI_Recv(Vec<Vec<TN> > &vec,int rank,int tag=201) {
    MPI_Status status;
    Vec<TN> a;
    MPI_Probe(rank,MPI_ANY_TAG,MPI_COMM_WORLD, &status );
//    a.resize( status.st_length / sizeof(TN));
    a.resize( status._count / sizeof(TN));
    MPI_Recv( a.ptr(), a.size(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &status );
    int size=0;
    for( int i=0 ;i<vec.size() ;i++ ){
        vec[i]=a[range(vec[i].size())+size ];
        size+=vec[i].size();
    }
}


/** \ingroup MPI_Recv
\brief Recevoir une matrice full de la classe LMT - support du tag

Il suffit de lancer la commande MPI_Recv(matrice, rank) pour recevoir la matrice envoy�e dans matrice par la machine rank. Il est n�cessaire de connaitre le type � recevoir ainsi que la taille.
*/
template<class TN>
void MPI_Recv(Mat<TN> &mat,int rank,int tag=201) {
    MPI_Status status;
    MPI_Recv( mat.data.ptr(), mat.nb_cols()*mat.nb_rows(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &status );
}


/** \ingroup MPI_Recv
\brief Recevoir une matrice sparse de la classe LMT - support du tag

Il suffit de lancer la commande MPI_Recv(matrice, rank) pour recevoir la matrice envoy�e dans matrice par la machine rank. Il est n�cessaire de connaitre le type � recevoir. La matrice sera resiz� automatiquement.
*/
template<class TN, class TO>
void MPI_Recv(Mat<TN,TO,SparseLine<Col> > &mat,int rank,int tag=201) {
    MPI_Status status;
    Vec<TN> a;
    MPI_Probe(rank,MPI_ANY_TAG,MPI_COMM_WORLD, &status );
//    a.resize( status.st_length / sizeof(TN));
    a.resize( status._count / sizeof(TN));
    MPI_Recv( a.ptr(), a.size(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &status );
    mat.resize((int) a[0], (int) a[1]);int pos=2;
    for( int i=0 ;i<mat.data.size() ;i++ ){
         mat.data[i].indices.resize((int) a[pos]);
         mat.data[i].data.resize((int) a[pos]);pos+=1;
         for( int j=0 ; j<mat.data[i].indices.size();j++ ){
             mat.data[i].indices[j]=(unsigned int) a[pos];pos+=1;
         }
          for( int j=0 ; j<mat.data[i].indices.size();j++ ){
             mat.data[i].data[j]=a[pos];pos+=1;
         }
    }
}


/*/** \ingroup MPI_Recv
\brief Recevoir un vecteur de la classe Time

Il suffit de lancer la commande MPI_Recv(vecteur, rank) pour recevoir le vecteur de Time envoy� dans vecteur par la machine rank. Il est n�cessaire de connaitre le type � recevoir ainsi que le pattern du vecteur r�ception (c'est-�-dire la taille du vecteur t et la taille des sous-vecteurs).
*/
// template<class TN>
// void MPI_Recv(Vec< Time<TN> > &t,int rank, Vec<string> toget){
// MPI_Recv(t,rank,toget,201);
// }
/*/** \ingroup MPI_Recv
\brief Recevoir un vecteur de la classe Time - support du tag
*/
// template<class TN>
// void MPI_Recv(Vec< Time<TN> > &t,int rank, Vec<string> toget,int tag){
//     MPI_Status status;
//     Vec<TN> a;
//     MPI_Probe(rank,MPI_ANY_TAG,MPI_COMM_WORLD, &status );
//     a.resize( t[0].W1.size() * toget.size() * t.size());
//     MPI_Recv( a.ptr(), a.size(), Mpi_type< typename Vec<TN>::template SubType<0> ::T >::res(), rank, tag, MPI_COMM_WORLD, &status );
//     int size=0;
//     for( int j =0 ;j < t.size() ; j++){
//         for( int i=0; i < toget.size() ;i++ ){
//             if (toget[i] == "Wd") {
//                 t[j].Wd=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Qd") {
//                 t[j].Qd=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wtemp") {
//                 t[j].Wtemp=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "F") {
//                 t[j].F=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wp") {
//                 t[j].Wp=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "W") {
//                 t[j].W=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Fchap") {
//                 t[j].Fchap=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wpchap") {
//                 t[j].Wpchap=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wchap") {
//                 t[j].Wchap=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "WtildeM") {
//                 t[j].WtildeM=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "oldF") {
//                 t[j].oldF=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "oldWp") {
//                 t[j].oldWp=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "oldW") {
//                 t[j].oldW=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "W1") {
//                 t[j].W1=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "W2") {
//                 t[j].W2=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "F1") {
//                 t[j].F1=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "F2") {
//                 t[j].F2=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wp1") {
//                 t[j].Wp1=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//             if (toget[i] == "Wp2") {
//                 t[j].Wp2=a[range(t[0].W1.size())+size];size+=t[0].W1.size();
//             }
//         }
//     }
// }



 /** \ingroup MPI
\defgroup MPI_Bcast
\brief Fonctions mettre une donnee sur tout les PRO

  */
/** \ingroup MPI_Bcast
\brief Bcast d un string

 */
template <class T>
void MPI_Bcast(T &str,int source){
    int size;
    size = str.size() + 1;
    MPI_Bcast(&size,1,MPI_INT,source,MPI_COMM_WORLD);
    char * buffer = new char[ size ];
    strncpy( buffer, str.c_str(), size );
    MPI_Bcast(buffer,size,MPI_CHAR,source,MPI_COMM_WORLD);
    str.assign(buffer);
}


#endif  //MPI_LMT_H
