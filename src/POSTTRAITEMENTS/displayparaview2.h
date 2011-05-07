//
// C++ Interface: displayparaview
//
// Description: 
//
//
// Author: Hugo LECLERC <leclerc@lmt.ens-cachan.fr>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef LMTDISPLAYPARAVIEW_H
#define LMTDISPLAYPARAVIEW_H

#include "mesh/write_mesh_vtk.h"
#include "implicitgeometry/shapedeclaration.h"
#include "containers/algo.h"

#include <string>
#include <fstream>
#include <sstream>

namespace LMT {

/**
open paraview and send data

@author Hugo LECLERC
*/
class DisplayParaview {
public:
    enum TypeField {Nodal,Elementary};
    DisplayParaview() {
        background_color = Vec<double,3>(0.2,0.2,0.3);
        init_xminmax = false;
    }
    ~DisplayParaview() {
    }
    
    template<class TM> void add_mesh(const TM &m,const std::string &prefix="tmp/paraview",const Vec<std::string> &display_fields=Vec<std::string>("all")) {
        std::ostringstream ss;
        ss << prefix << pvu_files.size() << ".vtu";
        std::string pvu_name( ss.str() );
        pvu_files.push_back(pvu_name);
        std::ofstream f(pvu_name.c_str());
        
        write_mesh_vtk<true>(f,m,display_fields);

        typename TM::Pvec xmi,xma;
        get_min_max( m.node_list, ExtractDM<pos_DM>(), xmi, xma );
        if ( m.node_list.size() )
            app_xminmax(xmi,xma);
    }
    template<class TS> void add_shape(const Shape<2,TS> &shape,unsigned grid_size,const std::string &prefix="tmp/paraview") {
        typedef typename Shape<2,TS>::Pvec Pvec;
        typedef typename Shape<2,TS>::TPen TPen;
        Pvec tdim_min,tdim_max;
        shape.get_dim_min_max(tdim_min,tdim_max);
        tdim_min -= 2.0/grid_size*(tdim_max-tdim_min);
        tdim_max += 2.0/grid_size*(tdim_max-tdim_min);
        app_xminmax(tdim_min,tdim_max);

        Vec<Pvec> points(grid_size*grid_size);
        Vec<TPen> pem(grid_size*grid_size);
        Pvec spacing = (tdim_max-tdim_min)/(grid_size-1);
        for (unsigned i=0,offset=0;i<grid_size;++i)
            for (unsigned j=0;j<grid_size;++j,++offset)
                points[offset] = tdim_min + Pvec(j,i) * spacing;
        shape.get_penetration( points.begin(), pem.begin(), points.size() );
        
        std::ostringstream ss;
        ss << prefix << vti_files.size() << ".vti";
        std::string vti_name( ss.str() );
        vti_files.push_back(vti_name);
        std::ofstream f(vti_name.c_str());
        
        f << "<?xml version='1.0'?>" << std::endl;
        f << "<VTKFile type='ImageData' version='0.1' byte_order='LittleEndian'>" << std::endl;
        f << "  <ImageData WholeExtent='0 " << grid_size-1 << " 0 " << grid_size-1
          << " 0 0' Origin='" << tdim_min << " 0' Spacing='" << spacing << " 0.01'>" << std::endl;
        f << "    <Piece Extent='0 " << grid_size-1 << " 0 " << grid_size-1 << " 0 0'>" << std::endl;
        f << "      <PointData Scalars='offset'>" << std::endl;
        f << "        <DataArray type='Float64' Name='offset' format='ascii' NumberOfComponents='1'>" << std::endl;
        for(unsigned i=0;i<pem.size();++i)
            f << pem[i].dist << ' ';
        //f.write((char *)scalars.begin(),scalars.size()*sizeof(float));
        f << std::endl;
        f << "        </DataArray>" << std::endl;
        f << "        <DataArray type='Float64' Name='normals' format='ascii' NumberOfComponents='3'>" << std::endl;
        for(unsigned i=0;i<pem.size();++i)
            f << pem[i].normal[0] << ' ' << pem[i].normal[1] << " 0 ";
        //f.write((char *)scalars.begin(),scalars.size()*sizeof(float));
        f << std::endl;
        f << "        </DataArray>" << std::endl;
        f << "      </PointData>" << std::endl;
        f << "      <CellData>" << std::endl;
        f << "      </CellData>" << std::endl;
        f << "    </Piece>" << std::endl;
        f << "  </ImageData>" << std::endl;
        f << "</VTKFile>" << std::endl;
    }
    void set_field_to_display(const std::string &name,TypeField type) {
        field_to_display = name;
        type_field_to_display = type;
    }
    void exec(bool all_mesh=true,const std::string &prefix="tmp/paraview") {
        std::string tmp_file = prefix + ".pvs";
        std::ofstream pvs( tmp_file.c_str() );
        pvs << "" << std::endl;
        pvs << "# ParaView State Version 1.8" << std::endl;
        pvs << "" << std::endl;
        pvs << "set kw(vtkApplication) [$Application GetMainWindow]" << std::endl;
        pvs << "set kw(vtkMainWin) [$kw(vtkApplication) GetMainView]" << std::endl;
        // pvs << "set kw(vtkTemp322) [$kw(vtkApplication) GetAnimationInterface]" << std::endl;
        //pvs << "[$kw(vtkApplication) GetRotateCameraButton] SetState 1" << std::endl;
        //pvs << "$kw(vtkApplication) ChangeInteractorStyle 1" << std::endl;
        for(unsigned i=0;i<( all_mesh ? pvu_files.size() : min(pvu_files.size(),(unsigned)1) );++i) {
            pvs << "set kw(vtkgr" << i << ") [$kw(vtkApplication) InitializeReadCustom \"XMLUnstructuredGridReader\" \""
                << pvu_files[i] << "\"]" << std::endl;
            pvs << "$kw(vtkApplication) ReadFileInformation $kw(vtkgr" << i << ") \"" << pvu_files[i] << "\"" << std::endl;
            pvs << "$kw(vtkApplication) FinalizeRead $kw(vtkgr" << i << ") \"" << pvu_files[i] << "\"" << std::endl;
            pvs << "set kw(vtkTemp534) [$kw(vtkgr" << i << ") GetPVWidget {Filename}]" << std::endl;
            pvs << "$kw(vtkTemp534) SetValue \"" << pvu_files[i] << "\"" << std::endl;
            pvs << "$kw(vtkgr" << i << ") AcceptCallback" << std::endl;
            if ( field_to_display.size() ) {
                if ( type_field_to_display == Nodal )
                    pvs << "[$kw(vtkgr"<<i<<") GetPVOutput] ColorByPointField {" << field_to_display << "} 1" << std::endl;
                else
                    pvs << "[$kw(vtkgr"<<i<<") GetPVOutput] ColorByCellField {" << field_to_display << "} 1" << std::endl;
            }
            pvs << "" << std::endl;
        }
        for(unsigned i=0;i<vti_files.size();++i) {
            pvs << "set kw(vtkTemp525) [$kw(vtkApplication) InitializeReadCustom \"XMLImageDataReader\" \""
                << vti_files[i] << "\"]" << std::endl;
            pvs << "$kw(vtkApplication) ReadFileInformation $kw(vtkTemp525) \"" << vti_files[i] << "\"" << std::endl;
            pvs << "$kw(vtkApplication) FinalizeRead $kw(vtkTemp525) \"" << vti_files[i] << "\"" << std::endl;
            pvs << "set kw(vtkTemp534) [$kw(vtkTemp525) GetPVWidget {Filename}]" << std::endl;
            pvs << "$kw(vtkTemp534) SetValue \"" << vti_files[i] << "\"" << std::endl;
            pvs << "set kw(vtkTemp529) [$kw(vtkTemp525) GetPVWidget {PointSelection}]" << std::endl;
            pvs << "$kw(vtkTemp529) SetArrayStatus {offset} 1" << std::endl;
            pvs << "set kw(vtkTemp533) [$kw(vtkTemp525) GetPVWidget {CellSelection}]" << std::endl;
            pvs << "$kw(vtkTemp525) AcceptCallback" << std::endl;
            pvs << "set pvDisp(vtkTemp525) [$kw(vtkTemp525) GetPartDisplay]" << std::endl;
            pvs << "$pvDisp(vtkTemp525) SetOpacity 0.3" << std::endl;
            pvs << "[$kw(vtkTemp525) GetPVOutput] ColorByPointField {offset} 1" << std::endl;
            pvs << "set kw(vtkTemp602) [$kw(vtkApplication) CreatePVSource Contour]" << std::endl;
            pvs << "set kw(vtkTemp608) [$kw(vtkTemp602) GetPVWidget {Contour Values}]" << std::endl;
            pvs << "$kw(vtkTemp608) AddValue 0" << std::endl;
            pvs << "$kw(vtkTemp602) AcceptCallback" << std::endl;
        }
        if ( field_to_display.size() ) {
            pvs << "set kw(vtkTemp571) [$kw(vtkApplication) GetPVColorMap {" << field_to_display << "} 1]" << std::endl;
            // pvs << "$kw(vtkTemp571) SetScalarRange 0 11" << std::endl;
            pvs << "$kw(vtkTemp571) SetStartHSV 0.6667 1 1" << std::endl;
            pvs << "$kw(vtkTemp571) SetEndHSV 0 1 1" << std::endl;
            pvs << "$kw(vtkTemp571) VectorModeMagnitudeCallback" << std::endl;
            pvs << "$kw(vtkTemp571) SetScalarBarVisibility 1" << std::endl;
        }
        pvs << "$kw(vtkMainWin) SetRendererBackgroundColor "
            << background_color[0] << " " 
            << background_color[1] << " " 
            << background_color[2] << std::endl;
        pvs << "$kw(vtkMainWin) ParallelProjectionOff" << std::endl;
    
        pvs << "$kw(vtkApplication) SetCenterOfRotation "
            << (xmin[0]+xmax[0])/2.0 << " " 
            << (xmin[1]+xmax[1])/2.0 << " " 
            << (xmin[2]+xmax[2])/2.0 << std::endl;
        pvs << "$kw(vtkApplication) ResetCameraCallback" << std::endl;
        //    pvs << "$kw(vtkMainWin) Render" << std::endl;
        pvs.close();
        
        system( ("paraview "+tmp_file).c_str() );
    }
private:
    template<class PV> void app_xminmax(const PV &xmi,const PV &xma) {
        if ( init_xminmax ) {
            for(unsigned i=0;i<min(xmi.size(),(unsigned)3);++i) {
                xmin[i] = min(xmin[i],xmi[i]); xmax[i] = max(xmax[i],xma[i]);
            }
        }
        else {
            xmin = 0.0; xmax = 0.0;
            for(unsigned i=0;i<min(xmi.size(),(unsigned)3);++i) {
                xmin[i] = xmi[i]; xmax[i] = xma[i];
            }
            init_xminmax = true;
        }
    }

    Vec<std::string> pvu_files;
    Vec<std::string> vti_files;
    std::string field_to_display;
    TypeField type_field_to_display;
    Vec<double,3> background_color;
    Vec<double,3> xmin,xmax;
    bool init_xminmax;
};

/**
 * declare a DisplayParaview instance, push m and interact with user
 * @param m 
 * @param nodal_field_to_display 
 */
template<class TM> void display_mesh(const TM &m,const char *nodal_field_to_display="",DisplayParaview::TypeField type_field_to_display=DisplayParaview::Nodal) {
    DisplayParaview dp;
    dp.add_mesh(m);
    dp.set_field_to_display( nodal_field_to_display, type_field_to_display );
    dp.exec();
}

/**
 * usefull to get several Windows (using apply_mt)
 */
struct DpExec { void operator()(DisplayParaview &dp,unsigned i) const { dp.exec(false,"tmp/conf"+to_string(i)); } };

};

#endif
