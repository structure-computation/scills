#ifndef CREATE_FILE_PVD
#define CREATE_FILE_PVD

#include "../DEFINITIONS/SavingData.h"
#include "../DEFINITIONS/TimeData.h"
#include "../DEFINITIONS/LatinData.h"
#include "../DEFINITIONS/MultiResolutionData.h"
#include "../DEFINITIONS/Process.h"

#include <fstream>


inline void create_file_pvtu(Process &process, const Sc2String prefix="sst_", unsigned pt=0, unsigned pt1=0){
    Vec<Sc2String> display_fields;
    unsigned typepvtu=1; //defaut quand on sort les inter ou sst et contact
    if (prefix=="sst_" or prefix=="inter_" or prefix=="contact_") typepvtu=0;
    
    Sc2String name_data = prefix+process.affichage->name_data;
    std::ostringstream ss2;
    ss2<<process.affichage->repertoire_save<<prefix<< ( (typepvtu==0) ? process.affichage->name_data : "" ) << ( (typepvtu==0) ? "_" : "" ) << ( (pt1!=0) ?  to_string(pt1) :"" ) << ( (pt1!=0) ? "_" :"" ) << ( (typepvtu==0) ?  to_string(pt) :"" ) << ".pvtu";
    
    Sc2String nom(ss2.str());
    std::ofstream os(nom.c_str());
    os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl;
    os << "<PUnstructuredGrid GhostLevel=\"0\">"<<std::endl;
    //ajout des donnees aux noeuds
    os << "\t<PPointData Vectors=\"vit\">"<<std::endl;
    
    std::ostringstream ss1;
    ss1<<process.affichage->repertoire_save<<((typepvtu==0)? prefix+process.affichage->name_data:process.affichage->type_affichage)<<"_"<<((pt1!=0)?process.parallelisation->rank:1)<<"_"<<((pt1!=0)? to_string(pt1)+"_":"")<<((typepvtu==0)?"1":"0")<<".vtu";
    Sc2String nomvtu(ss1.str());
    ifstream is(nomvtu.c_str());
    unsigned nbdatapoint=1;
    unsigned nbdatacell=1;
    if (is) {
        XmlNode n;
        n.parse_file( nomvtu.c_str() );
        XmlNode n1=n.get_element("UnstructuredGrid");
        XmlNode n2=n1.get_element("Piece");
        XmlNode npoints = n2.get_element("PointData");
        nbdatapoint=npoints.nb_elements("DataArray");
        for(unsigned i=0;i<nbdatapoint;i++){
            XmlNode ndata=npoints.get_element("DataArray",i);
            os << "                <DataArray Name='" << ndata.get_attribute("Name") << "' NumberOfComponents='" << ndata.get_attribute("NumberOfComponents") << "' type='" << ndata.get_attribute("type")  << "'/>"<<std::endl;
        }
    }
    os << "\t</PPointData>"<<std::endl;
    //ajout des donnees aux cellules
    os << "\t<PCellData Vectors=\"vectors\">"<<std::endl;
    if (is){
        XmlNode n;
        n.parse_file( nomvtu.c_str() );
        XmlNode n1=n.get_element("UnstructuredGrid");
        XmlNode n2=n1.get_element("Piece");
        XmlNode ncell = n2.get_element("CellData");
        nbdatacell=ncell.nb_elements("DataArray");
        for(unsigned i=0;i<nbdatacell;i++){
            XmlNode ndata=ncell.get_element("DataArray",i);
            os << "                <DataArray Name='" << ndata.get_attribute("Name") << "' NumberOfComponents='" << ndata.get_attribute("NumberOfComponents") << "' type='" << ndata.get_attribute("type")  << "'/>"<<std::endl;
        }
    }
    os << "\t</PCellData>"<<std::endl;
    //type de point
    os << "\t<PPoints>"<<std::endl;
    os << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\""<<DIM<<"\"/>"<<std::endl;
    os << "\t</PPoints>"<<std::endl;
    
    //ajout des parties
    for(unsigned i=1;i<(unsigned)process.parallelisation->size;i++){
        os<<"<Piece Source=\""<<((typepvtu==0)? prefix+process.affichage->name_data:process.affichage->type_affichage)<<"_"<<i<<"_"<<((pt1!=0)? to_string(pt1)+"_":"")<<((typepvtu==0)?to_string(pt):"0")<<".vtu\"/>"<<std::endl;
    }
    
    os << "</PUnstructuredGrid>\n</VTKFile>"<<std::endl;
    if ( nbdatapoint == 0 and nbdatacell == 0) int tmp=system(("rm "+nom ).c_str());
}


inline void create_file_pvd(Process &process,const Sc2String prefix="sst_",unsigned pt=0){
   unsigned nbfiles = process.temps->nbpastemps;
   Sc2String name_data = process.affichage->repertoire_save+prefix+process.affichage->name_data;
   if (prefix=="contact_") nbfiles=process.latin->iter+1;
   //nom du fichier pvd
   std::ostringstream spvd;
   if (prefix=="contact_") spvd<<process.affichage->repertoire_save<<prefix<<process.affichage->name_data<<"_"<<pt<<".pvd";
   else spvd<<process.affichage->repertoire_save<<prefix<<process.affichage->name_data<<".pvd";
   Sc2String namepvd( spvd.str() );
   std::ofstream file_pvd;
   file_pvd.open( (namepvd).c_str(), ofstream::out);
   
   //ecriture entete
   file_pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;
   
   //ecriture des lignes correspondant aux fichier a lire
   for(unsigned i=1;i<=nbfiles;++i){
        if (process.parallelisation->is_multi_cpu() and prefix!="contact_") create_file_pvtu(process,prefix,i);
        std::ostringstream ss;
        if (pt == 0) {
            if (process.parallelisation->is_multi_cpu())
                ss << prefix << process.affichage->name_data << "_" << i <<  ".pvtu";
            else
                ss << prefix << process.affichage->name_data << "_" << i <<  ".vtu";
        } else {
            if (process.parallelisation->is_multi_cpu())
                ss << prefix << process.affichage->name_data << "_" << pt << "_" << i-1 <<  ".pvtu";
            else
                ss << prefix << process.affichage->name_data << "_" << pt << "_" << i-1 <<  ".vtu";
        }
        Sc2String nom1( ss.str() );
        std::ostringstream ss2;
        if (pt == 0)
            ss2 << " <DataSet timestep=\" " << i << "\" group=\"\" part=\"0\"\n\t\t file=\"" << nom1 << "\"/>\n" ;
        else
            ss2 << " <DataSet timestep=\" " << i-1 << "\" group=\"\" part=\"0\"\n\t\t file=\"" << nom1 << "\"/>\n" ;
        Sc2String ligne1(ss2.str());
        file_pvd << ligne1;
   }
   //fin du fichier pvd
   file_pvd << " </Collection> \n </VTKFile>";
}

//creation d'un fichier unique permettant d'acceder a tous les pas de temps et toutes les donnees par processeur pour une resolution
inline void create_file_pvd(Process &process, DataUser &data_user, const Sc2String prefix="sst_bulk"){
    Sc2String save_directory=process.affichage->repertoire_save+"results/"+prefix+"/";
    
    //creation du nom et du fichier pvd
    std::ostringstream spvd;
    Sc2String name_multiresolution="";
    if(process.multiresolution->nb_calculs>1)
        name_multiresolution<<"resolution_"<<process.multiresolution->m<<"_";
    
    spvd<<process.affichage->repertoire_save<<"results/"<<name_multiresolution.c_str()<<prefix<<".pvd";
    Sc2String namepvd( spvd.str() );
    std::ofstream file_pvd;
    file_pvd.open( (namepvd).c_str(), ofstream::out);
    
    //ecriture entete
    file_pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;
    
    unsigned ini_proc=1;
    if(not process.parallelisation->is_multi_cpu())
        ini_proc=0;
    else
        ini_proc=1;
    //ecriture des lignes correspondant aux fichier a lire
    //boucle sur le nb de processeur, les step et les piquets de temps
    process.temps->pt_cur=0;
    for(unsigned i_step=0;i_step<process.temps->nb_step;i_step++){
        process.temps->step_cur=i_step;
        for(unsigned i_pt = 0 ; i_pt < process.temps->time_step[i_step].nb_time_step; i_pt++){
            process.temps->time_step[i_step].pt_cur=i_pt;
            process.temps->pt_cur+=1;
            process.temps->t_cur=process.temps->time_step[i_step].t_ini+(i_pt+1)*process.temps->time_step[i_step].dt;
            for(unsigned i_proc=ini_proc;i_proc<process.parallelisation->size;++i_proc){
                //nom du fichier a lire pour le piquet de temps et le processeur considere
                std::ostringstream ss;
                ss << save_directory << name_multiresolution.c_str()<<process.affichage->name_data << "_proc_"<<i_proc<<"_time_"<<process.temps->pt_cur <<  ".vtu";
                Sc2String nom1( ss.str() );
                //designation du dataset sous paraview (valeur du piquet de temps, processeur)    
                std::ostringstream ss2;
                ss2 << " <DataSet timestep=\" " << process.temps->t_cur << "\" group=\"\" part=\" "<<i_proc<<" \"\n\t\t file=\"" << nom1 << "\"/>\n" ;
                Sc2String ligne1(ss2.str()); 
                file_pvd << ligne1;
            }          
        }
    }
    //fin du fichier pvd
    file_pvd << " </Collection> \n </VTKFile>";   
}

//creation d'un fichier unique permettant d'acceder a tous les pas de temps et toutes les donnees par processeur pour une resolution
inline void create_file_pvd_geometry(Process &process, DataUser &data_user, const Sc2String prefix="sst_bulk"){
   Sc2String save_directory=process.affichage->repertoire_save+"results/"+prefix+"/";
    
   //creation du nom et du fichier pvd
   std::ostringstream spvd;
   
   spvd<<process.affichage->repertoire_save<<"results/"<<prefix<<".pvd";
   Sc2String namepvd( spvd.str() );
   
   std::ofstream file_pvd;
   file_pvd.open( (namepvd).c_str(), std::ofstream::out);
   
   //ecriture entete
   file_pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<std::endl;
   
   unsigned ini_proc=1;
   if(not process.parallelisation->is_multi_cpu())
       ini_proc=0;
   else
       ini_proc=1;
   //ecriture des lignes correspondant aux fichier a lire
    for(unsigned i_proc=ini_proc;i_proc<process.parallelisation->size;++i_proc){
        //nom du fichier a lire pour le piquet de temps et le processeur considere
        std::ostringstream ss;
        ss << process.affichage->repertoire_save<<"results/"<<prefix<< "_proc_"<<i_proc <<  ".vtu";
        Sc2String nom1( ss.str() );
        //designation du dataset sous paraview (valeur du piquet de temps, processeur)    
        std::ostringstream ss2;
        ss2 << " <DataSet timestep=\"1\"" << " group=\"\" part=\" "<<i_proc<<" \"\n\t\t file=\"" << nom1 << "\"/>\n" ;
        Sc2String ligne1(ss2.str()); 
        file_pvd << ligne1;
    }          
   //fin du fichier pvd
   file_pvd << " </Collection> \n </VTKFile>";   
}

#endif //CREATE_FILE_PVD
