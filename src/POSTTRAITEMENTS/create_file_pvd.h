#ifndef CREATE_FILE_PVD
#define CREATE_FILE_PVD

#include "definition_PARAM_AFFICHAGE.h"
#include "definition_PARAM_TEMPS.h"
#include "definition_PARAM.h"
#include "definition_PARAM_LATIN.h"


inline void create_file_pvtu(Param &process, const string prefix="sst_", unsigned pt=0, unsigned pt1=0){
  int dim=process.dim;
  Vec<string> display_fields;
  unsigned typepvtu=1; //defaut quand on sort les inter ou sst et contact
  if (prefix=="sst_" or prefix=="inter_" or prefix=="contact_") typepvtu=0;
  
  string name_data = prefix+process.affichage->name_data;
  std::ostringstream ss2;
  ss2<<process.affichage->repertoire_save<<prefix<< ( (typepvtu==0) ? process.affichage->name_data : "" ) << ( (typepvtu==0) ? "_" : "" ) << ( (pt1!=0) ?  to_string(pt1) :"" ) << ( (pt1!=0) ? "_" :"" ) << ( (typepvtu==0) ?  to_string(pt) :"" ) << ".pvtu";
  
  string nom(ss2.str());
  std::ofstream os(nom.c_str());
  os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  os << "<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
  //ajout des donnees aux noeuds
  os << "\t<PPointData Vectors=\"depl\">"<<endl;
  
  std::ostringstream ss1;
  ss1<<process.affichage->repertoire_save<<((typepvtu==0)? prefix+process.affichage->name_data:process.affichage->type_affichage)<<"_"<<((pt1!=0)?process.rank:1)<<"_"<<((pt1!=0)? to_string(pt1)+"_":"")<<((typepvtu==0)?"1":"0")<<".vtu";
  string nomvtu(ss1.str());
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
    os << "                <DataArray Name='" << ndata.get_attribute("Name") << "' NumberOfComponents='" << ndata.get_attribute("NumberOfComponents") << "' type='" << ndata.get_attribute("type")  << "'/>"<<endl;
  }
  }
  os << "\t</PPointData>"<<endl;
  //ajout des donnees aux cellules
  os << "\t<PCellData Vectors=\"vectors\">"<<endl;
  if (is){
      XmlNode n;
      n.parse_file( nomvtu.c_str() );
      XmlNode n1=n.get_element("UnstructuredGrid");
      XmlNode n2=n1.get_element("Piece");
      XmlNode ncell = n2.get_element("CellData");
  nbdatacell=ncell.nb_elements("DataArray");
  for(unsigned i=0;i<nbdatacell;i++){
    XmlNode ndata=ncell.get_element("DataArray",i);
    os << "                <DataArray Name='" << ndata.get_attribute("Name") << "' NumberOfComponents='" << ndata.get_attribute("NumberOfComponents") << "' type='" << ndata.get_attribute("type")  << "'/>"<<endl;
  }
  }
  os << "\t</PCellData>"<<endl;
 //type de point
  os << "\t<PPoints>"<<endl;
  os << "\t\t<DataArray type=\"Float64\" NumberOfComponents=\""<<dim<<"\"/>"<<endl;
  os << "\t</PPoints>"<<endl;
 //ajout des parties
  for(unsigned i=1;i<(unsigned)process.size;i++){
    os<<"<Piece Source=\""<<((typepvtu==0)? prefix+process.affichage->name_data:process.affichage->type_affichage)<<"_"<<i<<"_"<<((pt1!=0)? to_string(pt1)+"_":"")<<((typepvtu==0)?to_string(pt):"0")<<".vtu\"/>"<<endl;
  }
  os << "</PUnstructuredGrid>\n</VTKFile>"<<endl;
  if ( nbdatapoint == 0 and nbdatacell == 0) system(("rm "+nom ).c_str());
}


inline void create_file_pvd(Param &process,const string prefix="sst_",unsigned pt=0){
   unsigned nbfiles = process.temps->nbpastemps;
   string name_data = process.affichage->repertoire_save+prefix+process.affichage->name_data;
   if (prefix=="contact_") nbfiles=process.latin->iter+1;
   //nom du fichier pvd
   std::ostringstream spvd;
   if (prefix=="contact_") spvd<<process.affichage->repertoire_save<<prefix<<process.affichage->name_data<<"_"<<pt<<".pvd";
   else spvd<<process.affichage->repertoire_save<<prefix<<process.affichage->name_data<<".pvd";
   string namepvd( spvd.str() );
   std::ofstream file_pvd;
   file_pvd.open( (namepvd).c_str(), ofstream::out);
   
   //ecriture entete
   file_pvd << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" \n\t byte_order=\"LittleEndian\" \n\t compressor=\"vtkZLibDataCompressor\" >\n <Collection> " <<endl;
   
   //ecriture des lignes correspondant aux fichier a lire
   for(unsigned i=1;i<=nbfiles;++i){
     if (process.size>1 and prefix!="contact_") create_file_pvtu(process,prefix,i);
      std::ostringstream ss;
      if (pt == 0) {
      if (process.size > 1 ) ss << prefix << process.affichage->name_data << "_" << i <<  ".pvtu";
      else ss << prefix << process.affichage->name_data << "_" << i <<  ".vtu";
      } else {
        if (process.size > 1 ) ss << prefix << process.affichage->name_data << "_" << pt << "_" << i-1 <<  ".pvtu";
        else ss << prefix << process.affichage->name_data << "_" << pt << "_" << i-1 <<  ".vtu";
      }
      string nom1( ss.str() );
      std::ostringstream ss2;
      if (pt == 0) ss2 << " <DataSet timestep=\" " << i << "\" group=\"\" part=\"0\"\n\t\t file=\"" << nom1 << "\"/>\n" ;
      else ss2 << " <DataSet timestep=\" " << i-1 << "\" group=\"\" part=\"0\"\n\t\t file=\"" << nom1 << "\"/>\n" ;
      string ligne1(ss2.str());
      file_pvd << ligne1;
   }
   //fin du fichier pvd
   file_pvd << " </Collection> \n </VTKFile>";
}



#endif //CREATE_FILE_PVD
