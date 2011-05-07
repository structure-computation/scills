using namespace LMT;
using namespace std;

inline void read_data_structure(Param &process, const XmlNode &n) {

   XmlNode nm = n.get_element("mesh");
   nm.get_attribute( "repertoire_des_maillages", process.structure->repertoire_des_maillages );
   nm.get_attribute( "nom_fichier_qualification_materiaux", process.structure->nom_fichier_qualification_materiaux );
   nm.get_attribute( "nom_des_maillages", process.structure->nom_des_maillages );
   nm.get_attribute( "nb_de_maillages", process.structure->nb_maillages );
   nm.get_attribute( "extension", process.structure->extension );
   //nm.get_attribute( "nom_des_maillages_CL", process.structure->nom_des_maillages_CL );
   //nm.get_attribute( "nb_de_maillages_CL", process.structure->nb_maillages_CL );
   //nm.get_attribute( "echelle", process.structure->scale );
   nm.get_attribute( "jeu_physique", process.structure->jeu_physique );
   if ( process.structure->jeu_physique == 1) {
      nm.get_attribute( "nom_des_maillages_jeu", process.structure->nom_maillages_jeu );
      string stringter = nm.get_attribute("inter_jeu");
      Vec<string> v=tokenize(stringter,';');
      process.structure->inter_jeu.resize(v.size());
      for(unsigned i=0;i<v.size();i++) {
         Vec<string> v2 = tokenize(v[i],' ');
         process.structure->inter_jeu[i].resize(2);
         process.structure->inter_jeu[i][0]=atoi(v2[0].c_str());
         process.structure->inter_jeu[i][1]=atoi(v2[1].c_str());
      }
   }

};
