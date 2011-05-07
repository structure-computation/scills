
//creation des fichiers de donnees contenant :
// 1 :  les pas de temps + iterations + numero du fichier correspondant + les numeros des interfaces degradees
// 2 : les contraintes et deformations moyennes
inline void create_file_output(std::ofstream &fptit, std::ofstream &f, AFFICHAGE &process){
   string file_pt_it="_num_it_pt";
   fptit.open ((process.repertoire_save+process.name_data+file_pt_it).c_str(), ofstream::out );
   fptit << "numfile   pastemps   iteration  num_interfaces_degradee"<<endl;
   
   string name_file1="_contraintes_deformation_moyennes";
   f.open( (process.repertoire_save+process.name_data+name_file1).c_str() , ofstream::out );
   f << "sigmoy[xx yy zz xy xz yz] epsmoy[xx yy zz xy xz yz]  volumetot"<<endl;
}

//stockage du resultat sous forme de fichier paraview
//modification du nom du fichier de sortie
inline void modif_name_paraviewfile(AFFICHAGE &process, unsigned numfile){
   std::ostringstream ss;
   ss << process.name_data << "_" << numfile <<  ".vtu";
   string name_file2( ss.str() );
   cout << "stockage resultat : " << name_file2 << endl;
   system( ("mv ./tmp/paraview0.vtu "+process.repertoire_save+name_file2).c_str() );
}

inline void modif_name_paraviewfile_inter(AFFICHAGE &process, unsigned numfile){
   std::ostringstream ss;
   ss << process.name_data << "_inter_" << numfile <<  ".vtu";
   string name_file2( ss.str() );
   cout << "stockage resultat : " << name_file2 << endl;
   system( ("mv ./tmp/paraview0.vtu "+process.repertoire_save+name_file2).c_str() );
}
