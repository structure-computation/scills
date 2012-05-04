#ifndef ALL_DECLARATIONS
#define ALL_DECLARATIONS

class Process;
class Sst;
class Interface;
class Boundary;

#include "COMPUTE/DataUser.h"
#include "../LMT/include/containers/vec.h"
#include "../LMT/include/containers/vecpointedvalues.h"
#include "../LMT/include/containers/mat.h"
using namespace LMT;
using namespace Metil;

typedef Vec<TYPEREEL,DIM> ScPoint;
typedef Vec<TYPEREEL,DIM*(DIM+1)/2> ScVoigtVector;
typedef Vec<Sst> ScSstVec;
typedef Vec<VecPointedValues<Sst> > ScSstRef;
typedef Vec<Interface> ScInterVec;
typedef Vec<VecPointedValues<Interface> > ScInterRef;
typedef Vec<Boundary> ScCLVec;


#endif