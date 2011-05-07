#ifndef LMT___SRC_MAILLAGE_MESHCARACINTER2_MESHCARAC_PY
#define LMT___SRC_MAILLAGE_MESHCARACINTER2_MESHCARAC_PY

#include "mesh/triangle.h"
#include "mesh/quad_8.h"
#include "mesh/bar_3.h"
#include "mesh/tetra.h"
#include "mesh/triangle_6.h"
#include "mesh/quad.h"
#include "mesh/nodalelement.h"
#include "mesh/bar.h"
namespace LMT {

template<unsigned nb_dim,unsigned nvi> class Meshcaracinter2;


#ifndef IFNDEF_pos_DM
#define IFNDEF_pos_DM
    struct pos_DM {};
#endif // IFNDEF_pos_DM

template<>
class Meshcaracinter2<1,0> { 
public:
    static const unsigned dim = 1;
    typedef double Tpos;
    typedef Vec<double,1> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<2,0> { 
public:
    static const unsigned dim = 2;
    typedef double Tpos;
    typedef Vec<double,2> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<3,0> { 
public:
    static const unsigned dim = 3;
    typedef double Tpos;
    typedef Vec<double,3> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<1,1> { 
public:
    static const unsigned dim = 1;
    typedef double Tpos;
    typedef Vec<double,1> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<2,1> { 
public:
    static const unsigned dim = 2;
    typedef double Tpos;
    typedef Vec<double,2> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<3,1> { 
public:
    static const unsigned dim = 3;
    typedef double Tpos;
    typedef Vec<double,3> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<1,2> { 
public:
    static const unsigned dim = 1;
    typedef double Tpos;
    typedef Vec<double,1> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<2,2> { 
public:
    static const unsigned dim = 2;
    typedef double Tpos;
    typedef Vec<double,2> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<3,2> { 
public:
    static const unsigned dim = 3;
    typedef double Tpos;
    typedef Vec<double,3> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<1,3> { 
public:
    static const unsigned dim = 1;
    typedef double Tpos;
    typedef Vec<double,1> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Tetra NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<2,3> { 
public:
    static const unsigned dim = 2;
    typedef double Tpos;
    typedef Vec<double,2> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Tetra NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};

template<>
class Meshcaracinter2<3,3> { 
public:
    static const unsigned dim = 3;
    typedef double Tpos;
    typedef Vec<double,3> Pvec;
    struct NodalStaticData {
        CARACDMEXTNAME(0,Pvec,pos,"m");
        static const unsigned nb_params = 1;
    };
    struct GlobalStaticData {
        VOIDDMSET;
    };
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0>
    struct ElementChoice {
        typedef void NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<0,0,0,inner> {
        typedef Tetra NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<3,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<2,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
    template<unsigned inner> struct ElementChoice<3,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
    };
};



};

#endif // LMT___SRC_MAILLAGE_MESHCARACINTER2_MESHCARAC_PY
