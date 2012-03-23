
//procedure permettant de stocker la geometrie des interfaces en indiquant les interfaces degradables testees
// sortie sous forme de fichier paraview : data_inter_testees_numfile_repI.vtu
// ou numfile est le numero du fichier et repI la configuration testee parmi les n possibles pour cette iteration
template<class TV2> void stockage_interfaces_testees(TV2 &Inter, AFFICHAGE &process, const unsigned numfile, Vec<unsigned> &repI,const int inddeg){
   //affich_INTER_testees(Inter,process.save,repI,inddeg);
   std::ostringstream ss;
   ss << process.name_data << "_inter_testee_" << numfile << ".vtu";
   Sc2String name_file( ss.str() );
   system( ("mv ./tmp/paraview0.vtu "+process.repertoire_save+name_file).c_str() );
}


