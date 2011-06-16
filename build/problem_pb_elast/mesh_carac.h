#ifndef Mesh_carac_pb_elast_HEADER
#define Mesh_carac_pb_elast_HEADER
#include "mesh/displayparaview.h"
#include "mesh/triangle.h"
#include "mesh/quad.h"
#include "mesh/tetra.h"
#include "mesh/wedge.h"
#include "mesh/hexa.h"
namespace LMT {

template<class TP,unsigned dim> struct Mesh_carac_pb_elast {};
#ifndef IFNDEF_ener_DM
#define IFNDEF_ener_DM
    struct ener_DM { static std::string name() { return "ener"; } };
#endif // IFNDEF_ener_DM

#ifndef IFNDEF_pos_DM
#define IFNDEF_pos_DM
    struct pos_DM { static std::string name() { return "pos"; } };
#endif // IFNDEF_pos_DM

#ifndef IFNDEF_numsst_DM
#define IFNDEF_numsst_DM
    struct numsst_DM { static std::string name() { return "numsst"; } };
#endif // IFNDEF_numsst_DM

#ifndef IFNDEF_f_vol_DM
#define IFNDEF_f_vol_DM
    struct f_vol_DM { static std::string name() { return "f_vol"; } };
#endif // IFNDEF_f_vol_DM

#ifndef IFNDEF_num_proc_skin_DM
#define IFNDEF_num_proc_skin_DM
    struct num_proc_skin_DM { static std::string name() { return "num_proc_skin"; } };
#endif // IFNDEF_num_proc_skin_DM

#ifndef IFNDEF_qtrans_DM
#define IFNDEF_qtrans_DM
    struct qtrans_DM { static std::string name() { return "qtrans"; } };
#endif // IFNDEF_qtrans_DM

#ifndef IFNDEF_deltaT_DM
#define IFNDEF_deltaT_DM
    struct deltaT_DM { static std::string name() { return "deltaT"; } };
#endif // IFNDEF_deltaT_DM

#ifndef IFNDEF_density_DM
#define IFNDEF_density_DM
    struct density_DM { static std::string name() { return "density"; } };
#endif // IFNDEF_density_DM

#ifndef IFNDEF_sigma_local_skin_DM
#define IFNDEF_sigma_local_skin_DM
    struct sigma_local_skin_DM { static std::string name() { return "sigma_local_skin"; } };
#endif // IFNDEF_sigma_local_skin_DM

#ifndef IFNDEF_typmat_skin_DM
#define IFNDEF_typmat_skin_DM
    struct typmat_skin_DM { static std::string name() { return "typmat_skin"; } };
#endif // IFNDEF_typmat_skin_DM

#ifndef IFNDEF_epsilon_DM
#define IFNDEF_epsilon_DM
    struct epsilon_DM { static std::string name() { return "epsilon"; } };
#endif // IFNDEF_epsilon_DM

#ifndef IFNDEF_num_proc_DM
#define IFNDEF_num_proc_DM
    struct num_proc_DM { static std::string name() { return "num_proc"; } };
#endif // IFNDEF_num_proc_DM

#ifndef IFNDEF_epsilon_skin_DM
#define IFNDEF_epsilon_skin_DM
    struct epsilon_skin_DM { static std::string name() { return "epsilon_skin"; } };
#endif // IFNDEF_epsilon_skin_DM

#ifndef IFNDEF_poisson_ratio_DM
#define IFNDEF_poisson_ratio_DM
    struct poisson_ratio_DM { static std::string name() { return "poisson_ratio"; } };
#endif // IFNDEF_poisson_ratio_DM

#ifndef IFNDEF_alpha_DM
#define IFNDEF_alpha_DM
    struct alpha_DM { static std::string name() { return "alpha"; } };
#endif // IFNDEF_alpha_DM

#ifndef IFNDEF_numsst_skin_DM
#define IFNDEF_numsst_skin_DM
    struct numsst_skin_DM { static std::string name() { return "numsst_skin"; } };
#endif // IFNDEF_numsst_skin_DM

#ifndef IFNDEF_elastic_modulus_DM
#define IFNDEF_elastic_modulus_DM
    struct elastic_modulus_DM { static std::string name() { return "elastic_modulus"; } };
#endif // IFNDEF_elastic_modulus_DM

#ifndef IFNDEF_typmat_DM
#define IFNDEF_typmat_DM
    struct typmat_DM { static std::string name() { return "typmat"; } };
#endif // IFNDEF_typmat_DM

#ifndef IFNDEF_f_vol_e_DM
#define IFNDEF_f_vol_e_DM
    struct f_vol_e_DM { static std::string name() { return "f_vol_e"; } };
#endif // IFNDEF_f_vol_e_DM

#ifndef IFNDEF_dep_DM
#define IFNDEF_dep_DM
    struct dep_DM { static std::string name() { return "dep"; } };
#endif // IFNDEF_dep_DM

#ifndef IFNDEF_resolution_DM
#define IFNDEF_resolution_DM
    struct resolution_DM { static std::string name() { return "resolution"; } };
#endif // IFNDEF_resolution_DM

#ifndef IFNDEF_sigma_skin_DM
#define IFNDEF_sigma_skin_DM
    struct sigma_skin_DM { static std::string name() { return "sigma_skin"; } };
#endif // IFNDEF_sigma_skin_DM

#ifndef IFNDEF_is_on_skin_DM
#define IFNDEF_is_on_skin_DM
    struct is_on_skin_DM { static std::string name() { return "is_on_skin"; } };
#endif // IFNDEF_is_on_skin_DM

#ifndef IFNDEF_sigma_DM
#define IFNDEF_sigma_DM
    struct sigma_DM { static std::string name() { return "sigma"; } };
#endif // IFNDEF_sigma_DM

#ifndef IFNDEF_sigma_von_mises_DM
#define IFNDEF_sigma_von_mises_DM
    struct sigma_von_mises_DM { static std::string name() { return "sigma_von_mises"; } };
#endif // IFNDEF_sigma_von_mises_DM

template<class TP>
struct Mesh_carac_pb_elast<TP,2> {
    typedef TP Tpos;
    static const unsigned dim = 2;
    typedef Vec<TP,2> Pvec;
    struct NodalStaticData {
        typedef Vec<double,2> T1;
        typedef Tpos T2;
        typedef Vec<Tpos,2> T0;
        NodalStaticData():qtrans(0),dep(0.0),is_on_skin(0) {}
        CARACDMEXTNAME( 0, T0, pos, "m" );
        CARACDMEXTNAME( 1, T1, qtrans, "mm" );
        CARACDMEXTNAME( 2, T0, dep, "mm" );
        CARACDMEXTNAME( 3, T2, is_on_skin, "" );
        static const unsigned nb_params = 4;
        void dm_data_set_field( const std::string field_name, Tpos value ) {
            if ( field_name == "pos" ) { pos = value; return; }
            if ( field_name == "dep" ) { dep = value; return; }
            if ( field_name == "is_on_skin" ) { is_on_skin = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
            if ( field_name == "pos" ) { pos = value; return; }
            if ( field_name == "dep" ) { dep = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
            if ( field_name == "is_on_skin" ) { return is_on_skin; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Tpos(0);
        }
        Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,1>();
        }
        Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
            if ( field_name == "pos" ) { return pos; }
            if ( field_name == "dep" ) { return dep; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,2>();
        }
        Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,3>();
        }
        Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,4>();
        }
        Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,5>();
        }
        Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,6>();
        }
        Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<1,1> >();
        }
        Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<2,2> >();
        }
        Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<3,3> >();
        }
        Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<4,4> >();
        }
        Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<5,5> >();
        }
        Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<6,6> >();
        }
        Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<1> >();
        }
        Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<2> >();
        }
        Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<3> >();
        }
    };
    struct GlobalStaticData {
        typedef Tpos T1;
        typedef Vec<Tpos,2> T0;
        GlobalStaticData():f_vol(0.0,0.0),deltaT(0),density(1),poisson_ratio(0.3),alpha(0),elastic_modulus(15e3),resolution(0) {}
        CARACDMEXTNAME( 0, T0, f_vol, "N/m^3" );
        CARACDMEXTNAME( 1, T1, deltaT, "degC" );
        CARACDMEXTNAME( 2, T1, density, "kg/mm^3" );
        CARACDMEXTNAME( 3, T1, poisson_ratio, "1" );
        CARACDMEXTNAME( 4, T1, alpha, "" );
        CARACDMEXTNAME( 5, T1, elastic_modulus, "N/mm^2" );
        CARACDMEXTNAME( 6, T1, resolution, "" );
        static const unsigned nb_params = 7;
        void dm_data_set_field( const std::string field_name, Tpos value ) {
            if ( field_name == "f_vol" ) { f_vol = value; return; }
            if ( field_name == "deltaT" ) { deltaT = value; return; }
            if ( field_name == "density" ) { density = value; return; }
            if ( field_name == "poisson_ratio" ) { poisson_ratio = value; return; }
            if ( field_name == "alpha" ) { alpha = value; return; }
            if ( field_name == "elastic_modulus" ) { elastic_modulus = value; return; }
            if ( field_name == "resolution" ) { resolution = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
            if ( field_name == "f_vol" ) { f_vol = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
            if ( field_name == "deltaT" ) { return deltaT; }
            if ( field_name == "density" ) { return density; }
            if ( field_name == "poisson_ratio" ) { return poisson_ratio; }
            if ( field_name == "alpha" ) { return alpha; }
            if ( field_name == "elastic_modulus" ) { return elastic_modulus; }
            if ( field_name == "resolution" ) { return resolution; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Tpos(0);
        }
        Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,1>();
        }
        Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
            if ( field_name == "f_vol" ) { return f_vol; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,2>();
        }
        Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,3>();
        }
        Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,4>();
        }
        Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,5>();
        }
        Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,6>();
        }
        Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<1,1> >();
        }
        Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<2,2> >();
        }
        Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<3,3> >();
        }
        Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<4,4> >();
        }
        Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<5,5> >();
        }
        Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<6,6> >();
        }
        Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<1> >();
        }
        Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<2> >();
        }
        Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<3> >();
        }
    };
    typedef Node<dim,Tpos,NodalStaticData> TNode;
    typedef ElementAncestor<TNode> EA;
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0> struct ElementChoice { typedef void NE; typedef DefaultBehavior BE; typedef VoidDMSet TData; };
    template<unsigned skin,unsigned inner> struct ElementChoice<0,skin,0,inner> { typedef Triangle NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Vec<Tpos,3>,1> T1;
            typedef Tpos T0;
            typedef Vec<Tpos,2> T2;
            TData():ener(0),numsst(0),epsilon(0),num_proc(0),typmat(0),f_vol_e(0.0,0.0),sigma(0),sigma_von_mises(0) {}
            CARACDMEXTNAME( 0, T0, ener, "N*mm" );
            CARACDMEXTNAME( 1, T0, numsst, "" );
            CARACDMEXTNAME( 2, T1, epsilon, "" );
            CARACDMEXTNAME( 3, T0, num_proc, "" );
            CARACDMEXTNAME( 4, T0, typmat, "" );
            CARACDMEXTNAME( 5, T2, f_vol_e, "N/m^3" );
            CARACDMEXTNAME( 6, T1, sigma, "N/mm^2" );
            CARACDMEXTNAME( 7, T0, sigma_von_mises, "N/mm^2" );
            static const unsigned nb_params = 8;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "ener" ) { ener = value; return; }
                if ( field_name == "numsst" ) { numsst = value; return; }
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "num_proc" ) { num_proc = value; return; }
                if ( field_name == "typmat" ) { typmat = value; return; }
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                if ( field_name == "sigma_von_mises" ) { sigma_von_mises = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "ener" ) { return ener; }
                if ( field_name == "numsst" ) { return numsst; }
                if ( field_name == "num_proc" ) { return num_proc; }
                if ( field_name == "typmat" ) { return typmat; }
                if ( field_name == "sigma_von_mises" ) { return sigma_von_mises; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                if ( field_name == "f_vol_e" ) { return f_vol_e; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "epsilon" ) { return epsilon[0]; } // hum
                if ( field_name == "sigma" ) { return sigma[0]; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<0,skin,1,inner> { typedef Quad NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Vec<Tpos,3>,1> T1;
            typedef Tpos T0;
            typedef Vec<Tpos,2> T2;
            TData():ener(0),numsst(0),epsilon(0),num_proc(0),typmat(0),f_vol_e(0.0,0.0),sigma(0),sigma_von_mises(0) {}
            CARACDMEXTNAME( 0, T0, ener, "N*mm" );
            CARACDMEXTNAME( 1, T0, numsst, "" );
            CARACDMEXTNAME( 2, T1, epsilon, "" );
            CARACDMEXTNAME( 3, T0, num_proc, "" );
            CARACDMEXTNAME( 4, T0, typmat, "" );
            CARACDMEXTNAME( 5, T2, f_vol_e, "N/m^3" );
            CARACDMEXTNAME( 6, T1, sigma, "N/mm^2" );
            CARACDMEXTNAME( 7, T0, sigma_von_mises, "N/mm^2" );
            static const unsigned nb_params = 8;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "ener" ) { ener = value; return; }
                if ( field_name == "numsst" ) { numsst = value; return; }
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "num_proc" ) { num_proc = value; return; }
                if ( field_name == "typmat" ) { typmat = value; return; }
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                if ( field_name == "sigma_von_mises" ) { sigma_von_mises = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "ener" ) { return ener; }
                if ( field_name == "numsst" ) { return numsst; }
                if ( field_name == "num_proc" ) { return num_proc; }
                if ( field_name == "typmat" ) { return typmat; }
                if ( field_name == "sigma_von_mises" ) { return sigma_von_mises; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                if ( field_name == "f_vol_e" ) { return f_vol_e; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "epsilon" ) { return epsilon[0]; } // hum
                if ( field_name == "sigma" ) { return sigma[0]; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<1,skin,0,inner> { typedef Bar NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,3> T1;
            typedef Tpos T0;
            TData():num_proc_skin(0),sigma_local_skin(0),typmat_skin(0),epsilon_skin(0),numsst_skin(0),sigma_skin(0) {}
            CARACDMEXTNAME( 0, T0, num_proc_skin, "" );
            CARACDMEXTNAME( 1, T1, sigma_local_skin, "" );
            CARACDMEXTNAME( 2, T0, typmat_skin, "" );
            CARACDMEXTNAME( 3, T1, epsilon_skin, "" );
            CARACDMEXTNAME( 4, T0, numsst_skin, "" );
            CARACDMEXTNAME( 5, T1, sigma_skin, "N/mm^2" );
            static const unsigned nb_params = 6;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "num_proc_skin" ) { num_proc_skin = value; return; }
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "typmat_skin" ) { typmat_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "numsst_skin" ) { numsst_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "num_proc_skin" ) { return num_proc_skin; }
                if ( field_name == "typmat_skin" ) { return typmat_skin; }
                if ( field_name == "numsst_skin" ) { return numsst_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "sigma_local_skin" ) { return sigma_local_skin; }
                if ( field_name == "epsilon_skin" ) { return epsilon_skin; }
                if ( field_name == "sigma_skin" ) { return sigma_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<2,skin,0,inner> { typedef NodalElement NE; typedef DefaultBehavior BE; typedef VoidDMSet TData;    };
};
#ifndef IFNDEF_ener_DM
#define IFNDEF_ener_DM
    struct ener_DM { static std::string name() { return "ener"; } };
#endif // IFNDEF_ener_DM

#ifndef IFNDEF_pos_DM
#define IFNDEF_pos_DM
    struct pos_DM { static std::string name() { return "pos"; } };
#endif // IFNDEF_pos_DM

#ifndef IFNDEF_numsst_DM
#define IFNDEF_numsst_DM
    struct numsst_DM { static std::string name() { return "numsst"; } };
#endif // IFNDEF_numsst_DM

#ifndef IFNDEF_f_vol_DM
#define IFNDEF_f_vol_DM
    struct f_vol_DM { static std::string name() { return "f_vol"; } };
#endif // IFNDEF_f_vol_DM

#ifndef IFNDEF_num_proc_skin_DM
#define IFNDEF_num_proc_skin_DM
    struct num_proc_skin_DM { static std::string name() { return "num_proc_skin"; } };
#endif // IFNDEF_num_proc_skin_DM

#ifndef IFNDEF_qtrans_DM
#define IFNDEF_qtrans_DM
    struct qtrans_DM { static std::string name() { return "qtrans"; } };
#endif // IFNDEF_qtrans_DM

#ifndef IFNDEF_deltaT_DM
#define IFNDEF_deltaT_DM
    struct deltaT_DM { static std::string name() { return "deltaT"; } };
#endif // IFNDEF_deltaT_DM

#ifndef IFNDEF_density_DM
#define IFNDEF_density_DM
    struct density_DM { static std::string name() { return "density"; } };
#endif // IFNDEF_density_DM

#ifndef IFNDEF_sigma_local_skin_DM
#define IFNDEF_sigma_local_skin_DM
    struct sigma_local_skin_DM { static std::string name() { return "sigma_local_skin"; } };
#endif // IFNDEF_sigma_local_skin_DM

#ifndef IFNDEF_typmat_skin_DM
#define IFNDEF_typmat_skin_DM
    struct typmat_skin_DM { static std::string name() { return "typmat_skin"; } };
#endif // IFNDEF_typmat_skin_DM

#ifndef IFNDEF_epsilon_DM
#define IFNDEF_epsilon_DM
    struct epsilon_DM { static std::string name() { return "epsilon"; } };
#endif // IFNDEF_epsilon_DM

#ifndef IFNDEF_num_proc_DM
#define IFNDEF_num_proc_DM
    struct num_proc_DM { static std::string name() { return "num_proc"; } };
#endif // IFNDEF_num_proc_DM

#ifndef IFNDEF_epsilon_skin_DM
#define IFNDEF_epsilon_skin_DM
    struct epsilon_skin_DM { static std::string name() { return "epsilon_skin"; } };
#endif // IFNDEF_epsilon_skin_DM

#ifndef IFNDEF_poisson_ratio_DM
#define IFNDEF_poisson_ratio_DM
    struct poisson_ratio_DM { static std::string name() { return "poisson_ratio"; } };
#endif // IFNDEF_poisson_ratio_DM

#ifndef IFNDEF_alpha_DM
#define IFNDEF_alpha_DM
    struct alpha_DM { static std::string name() { return "alpha"; } };
#endif // IFNDEF_alpha_DM

#ifndef IFNDEF_numsst_skin_DM
#define IFNDEF_numsst_skin_DM
    struct numsst_skin_DM { static std::string name() { return "numsst_skin"; } };
#endif // IFNDEF_numsst_skin_DM

#ifndef IFNDEF_elastic_modulus_DM
#define IFNDEF_elastic_modulus_DM
    struct elastic_modulus_DM { static std::string name() { return "elastic_modulus"; } };
#endif // IFNDEF_elastic_modulus_DM

#ifndef IFNDEF_typmat_DM
#define IFNDEF_typmat_DM
    struct typmat_DM { static std::string name() { return "typmat"; } };
#endif // IFNDEF_typmat_DM

#ifndef IFNDEF_f_vol_e_DM
#define IFNDEF_f_vol_e_DM
    struct f_vol_e_DM { static std::string name() { return "f_vol_e"; } };
#endif // IFNDEF_f_vol_e_DM

#ifndef IFNDEF_dep_DM
#define IFNDEF_dep_DM
    struct dep_DM { static std::string name() { return "dep"; } };
#endif // IFNDEF_dep_DM

#ifndef IFNDEF_resolution_DM
#define IFNDEF_resolution_DM
    struct resolution_DM { static std::string name() { return "resolution"; } };
#endif // IFNDEF_resolution_DM

#ifndef IFNDEF_sigma_skin_DM
#define IFNDEF_sigma_skin_DM
    struct sigma_skin_DM { static std::string name() { return "sigma_skin"; } };
#endif // IFNDEF_sigma_skin_DM

#ifndef IFNDEF_is_on_skin_DM
#define IFNDEF_is_on_skin_DM
    struct is_on_skin_DM { static std::string name() { return "is_on_skin"; } };
#endif // IFNDEF_is_on_skin_DM

#ifndef IFNDEF_sigma_DM
#define IFNDEF_sigma_DM
    struct sigma_DM { static std::string name() { return "sigma"; } };
#endif // IFNDEF_sigma_DM

#ifndef IFNDEF_sigma_von_mises_DM
#define IFNDEF_sigma_von_mises_DM
    struct sigma_von_mises_DM { static std::string name() { return "sigma_von_mises"; } };
#endif // IFNDEF_sigma_von_mises_DM

template<class TP>
struct Mesh_carac_pb_elast<TP,3> {
    typedef TP Tpos;
    static const unsigned dim = 3;
    typedef Vec<TP,3> Pvec;
    struct NodalStaticData {
        typedef Vec<Tpos,3> T0;
        typedef Vec<double,3> T1;
        typedef Tpos T2;
        NodalStaticData():qtrans(0),dep(0.0),is_on_skin(0) {}
        CARACDMEXTNAME( 0, T0, pos, "m" );
        CARACDMEXTNAME( 1, T1, qtrans, "mm" );
        CARACDMEXTNAME( 2, T0, dep, "mm" );
        CARACDMEXTNAME( 3, T2, is_on_skin, "" );
        static const unsigned nb_params = 4;
        void dm_data_set_field( const std::string field_name, Tpos value ) {
            if ( field_name == "pos" ) { pos = value; return; }
            if ( field_name == "dep" ) { dep = value; return; }
            if ( field_name == "is_on_skin" ) { is_on_skin = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
            if ( field_name == "pos" ) { pos = value; return; }
            if ( field_name == "dep" ) { dep = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
            if ( field_name == "is_on_skin" ) { return is_on_skin; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Tpos(0);
        }
        Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,1>();
        }
        Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,2>();
        }
        Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
            if ( field_name == "pos" ) { return pos; }
            if ( field_name == "dep" ) { return dep; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,3>();
        }
        Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,4>();
        }
        Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,5>();
        }
        Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,6>();
        }
        Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<1,1> >();
        }
        Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<2,2> >();
        }
        Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<3,3> >();
        }
        Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<4,4> >();
        }
        Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<5,5> >();
        }
        Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<6,6> >();
        }
        Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<1> >();
        }
        Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<2> >();
        }
        Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<3> >();
        }
    };
    struct GlobalStaticData {
        typedef Vec<Tpos,3> T0;
        typedef Tpos T1;
        GlobalStaticData():f_vol(0.0,0.0,0.0),deltaT(0),density(1),poisson_ratio(0.3),alpha(0),elastic_modulus(15e3),resolution(0) {}
        CARACDMEXTNAME( 0, T0, f_vol, "N/m^3" );
        CARACDMEXTNAME( 1, T1, deltaT, "degC" );
        CARACDMEXTNAME( 2, T1, density, "kg/mm^3" );
        CARACDMEXTNAME( 3, T1, poisson_ratio, "1" );
        CARACDMEXTNAME( 4, T1, alpha, "" );
        CARACDMEXTNAME( 5, T1, elastic_modulus, "N/mm^2" );
        CARACDMEXTNAME( 6, T1, resolution, "" );
        static const unsigned nb_params = 7;
        void dm_data_set_field( const std::string field_name, Tpos value ) {
            if ( field_name == "f_vol" ) { f_vol = value; return; }
            if ( field_name == "deltaT" ) { deltaT = value; return; }
            if ( field_name == "density" ) { density = value; return; }
            if ( field_name == "poisson_ratio" ) { poisson_ratio = value; return; }
            if ( field_name == "alpha" ) { alpha = value; return; }
            if ( field_name == "elastic_modulus" ) { elastic_modulus = value; return; }
            if ( field_name == "resolution" ) { resolution = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
            if ( field_name == "f_vol" ) { f_vol = value; return; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
        }
        Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
            if ( field_name == "deltaT" ) { return deltaT; }
            if ( field_name == "density" ) { return density; }
            if ( field_name == "poisson_ratio" ) { return poisson_ratio; }
            if ( field_name == "alpha" ) { return alpha; }
            if ( field_name == "elastic_modulus" ) { return elastic_modulus; }
            if ( field_name == "resolution" ) { return resolution; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Tpos(0);
        }
        Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,1>();
        }
        Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,2>();
        }
        Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
            if ( field_name == "f_vol" ) { return f_vol; }
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,3>();
        }
        Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,4>();
        }
        Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,5>();
        }
        Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Vec<Tpos,6>();
        }
        Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<1,1> >();
        }
        Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<2,2> >();
        }
        Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<3,3> >();
        }
        Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<4,4> >();
        }
        Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<5,5> >();
        }
        Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
            std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            return Mat<Tpos,Gen<6,6> >();
        }
        Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<1> >();
        }
        Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<2> >();
        }
        Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
            assert( 0 /*TODO*/ );
            return Mat<Tpos,Sym<3> >();
        }
    };
    typedef Node<dim,Tpos,NodalStaticData> TNode;
    typedef ElementAncestor<TNode> EA;
    template<unsigned nvi_to_subs,unsigned skin,unsigned num_sub_element,unsigned inner=0> struct ElementChoice { typedef void NE; typedef DefaultBehavior BE; typedef VoidDMSet TData; };
    template<unsigned skin,unsigned inner> struct ElementChoice<0,skin,0,inner> { typedef Tetra NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,3> T2;
            typedef Vec<Vec<Tpos,6>,1> T1;
            typedef Tpos T0;
            TData():ener(0),numsst(0),epsilon(0),num_proc(0),typmat(0),f_vol_e(0.0,0.0,0.0),sigma(0),sigma_von_mises(0) {}
            CARACDMEXTNAME( 0, T0, ener, "N*mm" );
            CARACDMEXTNAME( 1, T0, numsst, "" );
            CARACDMEXTNAME( 2, T1, epsilon, "" );
            CARACDMEXTNAME( 3, T0, num_proc, "" );
            CARACDMEXTNAME( 4, T0, typmat, "" );
            CARACDMEXTNAME( 5, T2, f_vol_e, "N/m^3" );
            CARACDMEXTNAME( 6, T1, sigma, "N/mm^2" );
            CARACDMEXTNAME( 7, T0, sigma_von_mises, "N/mm^2" );
            static const unsigned nb_params = 8;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "ener" ) { ener = value; return; }
                if ( field_name == "numsst" ) { numsst = value; return; }
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "num_proc" ) { num_proc = value; return; }
                if ( field_name == "typmat" ) { typmat = value; return; }
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                if ( field_name == "sigma_von_mises" ) { sigma_von_mises = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "ener" ) { return ener; }
                if ( field_name == "numsst" ) { return numsst; }
                if ( field_name == "num_proc" ) { return num_proc; }
                if ( field_name == "typmat" ) { return typmat; }
                if ( field_name == "sigma_von_mises" ) { return sigma_von_mises; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "f_vol_e" ) { return f_vol_e; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                if ( field_name == "epsilon" ) { return epsilon[0]; } // hum
                if ( field_name == "sigma" ) { return sigma[0]; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<0,skin,1,inner> { typedef Wedge NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,3> T2;
            typedef Vec<Vec<Tpos,6>,1> T1;
            typedef Tpos T0;
            TData():ener(0),numsst(0),epsilon(0),num_proc(0),typmat(0),f_vol_e(0.0,0.0,0.0),sigma(0),sigma_von_mises(0) {}
            CARACDMEXTNAME( 0, T0, ener, "N*mm" );
            CARACDMEXTNAME( 1, T0, numsst, "" );
            CARACDMEXTNAME( 2, T1, epsilon, "" );
            CARACDMEXTNAME( 3, T0, num_proc, "" );
            CARACDMEXTNAME( 4, T0, typmat, "" );
            CARACDMEXTNAME( 5, T2, f_vol_e, "N/m^3" );
            CARACDMEXTNAME( 6, T1, sigma, "N/mm^2" );
            CARACDMEXTNAME( 7, T0, sigma_von_mises, "N/mm^2" );
            static const unsigned nb_params = 8;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "ener" ) { ener = value; return; }
                if ( field_name == "numsst" ) { numsst = value; return; }
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "num_proc" ) { num_proc = value; return; }
                if ( field_name == "typmat" ) { typmat = value; return; }
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                if ( field_name == "sigma_von_mises" ) { sigma_von_mises = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "ener" ) { return ener; }
                if ( field_name == "numsst" ) { return numsst; }
                if ( field_name == "num_proc" ) { return num_proc; }
                if ( field_name == "typmat" ) { return typmat; }
                if ( field_name == "sigma_von_mises" ) { return sigma_von_mises; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "f_vol_e" ) { return f_vol_e; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                if ( field_name == "epsilon" ) { return epsilon[0]; } // hum
                if ( field_name == "sigma" ) { return sigma[0]; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<0,skin,2,inner> { typedef Hexa NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,3> T2;
            typedef Vec<Vec<Tpos,6>,1> T1;
            typedef Tpos T0;
            TData():ener(0),numsst(0),epsilon(0),num_proc(0),typmat(0),f_vol_e(0.0,0.0,0.0),sigma(0),sigma_von_mises(0) {}
            CARACDMEXTNAME( 0, T0, ener, "N*mm" );
            CARACDMEXTNAME( 1, T0, numsst, "" );
            CARACDMEXTNAME( 2, T1, epsilon, "" );
            CARACDMEXTNAME( 3, T0, num_proc, "" );
            CARACDMEXTNAME( 4, T0, typmat, "" );
            CARACDMEXTNAME( 5, T2, f_vol_e, "N/m^3" );
            CARACDMEXTNAME( 6, T1, sigma, "N/mm^2" );
            CARACDMEXTNAME( 7, T0, sigma_von_mises, "N/mm^2" );
            static const unsigned nb_params = 8;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "ener" ) { ener = value; return; }
                if ( field_name == "numsst" ) { numsst = value; return; }
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "num_proc" ) { num_proc = value; return; }
                if ( field_name == "typmat" ) { typmat = value; return; }
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                if ( field_name == "sigma_von_mises" ) { sigma_von_mises = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                if ( field_name == "f_vol_e" ) { f_vol_e = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                if ( field_name == "epsilon" ) { epsilon[0] = value; } // hum
                if ( field_name == "sigma" ) { sigma[0] = value; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "ener" ) { return ener; }
                if ( field_name == "numsst" ) { return numsst; }
                if ( field_name == "num_proc" ) { return num_proc; }
                if ( field_name == "typmat" ) { return typmat; }
                if ( field_name == "sigma_von_mises" ) { return sigma_von_mises; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                if ( field_name == "f_vol_e" ) { return f_vol_e; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                if ( field_name == "epsilon" ) { return epsilon[0]; } // hum
                if ( field_name == "sigma" ) { return sigma[0]; } // hum
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<1,skin,0,inner> { typedef Triangle NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,6> T1;
            typedef Tpos T0;
            TData():num_proc_skin(0),sigma_local_skin(0),typmat_skin(0),epsilon_skin(0),numsst_skin(0),sigma_skin(0) {}
            CARACDMEXTNAME( 0, T0, num_proc_skin, "" );
            CARACDMEXTNAME( 1, T1, sigma_local_skin, "" );
            CARACDMEXTNAME( 2, T0, typmat_skin, "" );
            CARACDMEXTNAME( 3, T1, epsilon_skin, "" );
            CARACDMEXTNAME( 4, T0, numsst_skin, "" );
            CARACDMEXTNAME( 5, T1, sigma_skin, "N/mm^2" );
            static const unsigned nb_params = 6;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "num_proc_skin" ) { num_proc_skin = value; return; }
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "typmat_skin" ) { typmat_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "numsst_skin" ) { numsst_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "num_proc_skin" ) { return num_proc_skin; }
                if ( field_name == "typmat_skin" ) { return typmat_skin; }
                if ( field_name == "numsst_skin" ) { return numsst_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                if ( field_name == "sigma_local_skin" ) { return sigma_local_skin; }
                if ( field_name == "epsilon_skin" ) { return epsilon_skin; }
                if ( field_name == "sigma_skin" ) { return sigma_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<1,skin,1,inner> { typedef Quad NE; typedef DefaultBehavior BE; 
        struct TData {
            typedef Vec<Tpos,6> T1;
            typedef Tpos T0;
            TData():num_proc_skin(0),sigma_local_skin(0),typmat_skin(0),epsilon_skin(0),numsst_skin(0),sigma_skin(0) {}
            CARACDMEXTNAME( 0, T0, num_proc_skin, "" );
            CARACDMEXTNAME( 1, T1, sigma_local_skin, "" );
            CARACDMEXTNAME( 2, T0, typmat_skin, "" );
            CARACDMEXTNAME( 3, T1, epsilon_skin, "" );
            CARACDMEXTNAME( 4, T0, numsst_skin, "" );
            CARACDMEXTNAME( 5, T1, sigma_skin, "N/mm^2" );
            static const unsigned nb_params = 6;
            void dm_data_set_field( const std::string field_name, Tpos value ) {
                if ( field_name == "num_proc_skin" ) { num_proc_skin = value; return; }
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "typmat_skin" ) { typmat_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "numsst_skin" ) { numsst_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,1> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,2> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,3> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,4> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,5> &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Vec<Tpos,6> &value ) {
                if ( field_name == "sigma_local_skin" ) { sigma_local_skin = value; return; }
                if ( field_name == "epsilon_skin" ) { epsilon_skin = value; return; }
                if ( field_name == "sigma_skin" ) { sigma_skin = value; return; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<1> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<2> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<3> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<4> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<5> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            void dm_data_set_field( const std::string field_name, const Mat<Tpos,Gen<6> > &value ) {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
            }
            Tpos dm_data_get_field( const std::string field_name, StructForType<Tpos> ) const {
                if ( field_name == "num_proc_skin" ) { return num_proc_skin; }
                if ( field_name == "typmat_skin" ) { return typmat_skin; }
                if ( field_name == "numsst_skin" ) { return numsst_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Tpos(0);
            }
            Vec<Tpos,1> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,1> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,1>();
            }
            Vec<Tpos,2> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,2> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,2>();
            }
            Vec<Tpos,3> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,3> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,3>();
            }
            Vec<Tpos,4> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,4> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,4>();
            }
            Vec<Tpos,5> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,5> > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,5>();
            }
            Vec<Tpos,6> dm_data_get_field( const std::string field_name, StructForType<Vec<Tpos,6> > ) const {
                if ( field_name == "sigma_local_skin" ) { return sigma_local_skin; }
                if ( field_name == "epsilon_skin" ) { return epsilon_skin; }
                if ( field_name == "sigma_skin" ) { return sigma_skin; }
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Vec<Tpos,6>();
            }
            Mat<Tpos,Gen<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<1,1> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<1,1> >();
            }
            Mat<Tpos,Gen<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<2,2> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<2,2> >();
            }
            Mat<Tpos,Gen<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<3,3> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<3,3> >();
            }
            Mat<Tpos,Gen<4> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<4,4> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<4,4> >();
            }
            Mat<Tpos,Gen<5> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<5,5> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<5,5> >();
            }
            Mat<Tpos,Gen<6> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Gen<6,6> > > ) const {
                std::cerr << "There is no variable named " << field_name << " in data struct" << std::endl;
                return Mat<Tpos,Gen<6,6> >();
            }
            Mat<Tpos,Sym<1> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<1> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<1> >();
            }
            Mat<Tpos,Sym<2> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<2> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<2> >();
            }
            Mat<Tpos,Sym<3> > dm_data_get_field( const std::string field_name, StructForType<Mat<Tpos,Sym<3> > > ) const {
                assert( 0 /*TODO*/ );
                return Mat<Tpos,Sym<3> >();
            }
        };
    };
    template<unsigned skin,unsigned inner> struct ElementChoice<2,skin,0,inner> { typedef Bar NE; typedef DefaultBehavior BE; typedef VoidDMSet TData;    };
    template<unsigned skin,unsigned inner> struct ElementChoice<3,skin,0,inner> { typedef NodalElement NE; typedef DefaultBehavior BE; typedef VoidDMSet TData;    };
};
} // namespace LMT
#endif // Mesh_carac_pb_elast_HEADER
