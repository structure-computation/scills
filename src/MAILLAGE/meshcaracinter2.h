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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<0,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
    };
    struct GlobalStaticData {
        VOIDDMSET;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<0, 0>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<1, 1>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<2, 2>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<3, 3>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<4, 4>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<5, 5>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Mat<double, LMT::Gen<6, 6>, LMT::Dense<LMT::Col> > &value, bool disp = true ) { assert(0); /*TODO*/ }
        template<class __G__> __G__ dm_data_get_field( const std::string field_name, StructForType<__G__>, bool disp = true ) const { assert( 0 /*TODO*/ ); return __G__( 0.0 );  }
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
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,0,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,0,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<3,0,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,0,inner> {
        typedef Triangle NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,1,inner> {
        typedef Quad NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,2,inner> {
        typedef Triangle_6 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<1,1,3,inner> {
        typedef Quad_8 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,0,inner> {
        typedef Bar NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<2,1,1,inner> {
        typedef Bar_3 NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
    template<unsigned inner> struct ElementChoice<3,1,0,inner> {
        typedef NodalElement NE;
        typedef DefaultBehavior BE;
        typedef VoidDMSet TData;
        void dm_data_set_field( const std::string field_name, Tpos value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value, bool disp = true ) { assert(0); /*TODO*/ }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value, bool disp = true ) { assert(0); /*TODO*/ }
    };
};



};

#endif // LMT___SRC_MAILLAGE_MESHCARACINTER2_MESHCARAC_PY
