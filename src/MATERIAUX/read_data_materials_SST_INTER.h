
#include "codegen/codegen.h"
#include "containers/basicops.h"
#include "DataUser.h"

using namespace Codegen;
using namespace LMT;
using namespace std;

/** \ingroup Materiaux
\brief Lecture des proprietes materiau par sst.
 
Le bloc <materials>  </materials> est lu dans le fichier xml.
 
On détermine ensuite le nombre de blocs <coefficients /> pour créer le vecteur de propriétés matériau (classe SstCarac).
 
L'"identificateur" permet de positionner le matériau lu dans le vecteur de propriétés. Le numéro affecté aux sous-structures correspond alors à cet identificateur.
 
Selon la dimension le champ "resolution" est lu ou non. Celui ci peut prendre les valeurs
\code resolution="deformation_plane" ou resolution="contrainte_plane" \endcode
 
Ensuite le champ "type" renseigne sur le comportement de la sous-structure. Selon la valeur obtenue, les différentes propriétés sont extraites du fichier xml (les champs nécessaires sont consultables dans les fichiers de formulations correspondants).
- \code type="isotrope" \endcode
- \code type="orthotrope" \endcode
- \code type="viscoelas" \endcode
- \code type="orthotrope_damage" \endcode
 
On aura donc les possibilités suivantes : 
-   On peut choisir un matériau isotrope  :
   \code 
   <coefficients type="isotrope" resolution="contrainte_plane" identificateur="0" name="fibre"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6"/>
   \endcode
-   ou bien un matériau orthotrope
   \code
   <coefficients type="orthotrope" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0"/>
   \endcode
-   ou encore un matériau orthotrope endommageable pour lequel on ajoute 
   \code
   <coefficients type="orthotrope_damage" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0" Yo="0.0961" Yop="0.0961" Ycp="10.8" Yc="10.8" b="2.5" Ysp="0.70"/>
   \endcode
-   ou encore un matériau viscoelastique (orthotrope ou isotrope)
   \code
   <coefficients type="isotrope" resolution="contrainte_plane" identificateur="0"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6" viscosite="0.1" />
   \endcode
*/




template<class TV3>
void read_material_properties(TV3 &matprop, Param &process, const DataUser &data_user) {

    unsigned nbmat = data_user.behaviour_materials.size();

    matprop.resize(nbmat);
    for(unsigned i=0;i<nbmat;++i) {
        matprop[i].id = data_user.behaviour_materials[i].id;
        matprop[i].type = data_user.behaviour_materials[i].type;
        matprop[i].comp = data_user.behaviour_materials[i].comp;
        if(data_user.dim == 2){
            if (data_user.behaviour_materials[i].resolution =="contrainte_plane")
                matprop[i].resolution=1;
            else if (data_user.behaviour_materials[i].resolution =="deformation_plane")
                matprop[i].resolution=0;
            else {
                std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
                assert(0);
            }
        }else{
            matprop[i].resolution=0;
        }

        std::vector<Ex> symbols;
        if (TV3::template SubType<0>::T::dim==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (TV3::template SubType<0>::T::dim==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }

        Vec<string> vstr;
        vstr.resize(data_user.dim);
        
        for(int d=0; d<data_user.dim; d++){
            vstr[d] = data_user.behaviour_bc_volume[2].step[0].CLv_step_prop[d] + " * " + data_user.behaviour_bc_volume[2].step[0].CLv_step_prop[6] ;
        }
        
        Vec<Ex> expr;
        expr.resize(TV3::template SubType<0>::T::dim);
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(vstr[d2],symbols);
        }
        Vec<double,TV3::template SubType<0>::T::dim> data;
        Ex::MapExNum var;
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= 0.;
        }
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            data[d2] = (double)expr[d2].subs_numerical(var);

        matprop[i].f_vol=data;
//         std::cout << "Pour le materiau  " << id << " : " << data << std::endl;

        Vec< TYPE > mat_prop_temp;
        mat_prop_temp.resize(data_user.behaviour_materials[i].mat_prop.size());
        for(int i_prop=0; i_prop<data_user.behaviour_materials[i].mat_prop.size(); i_prop++){
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_materials[i].mat_prop[i_prop],symbols);
            Ex::MapExNum var;
            for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= 0.;
            }
            mat_prop_temp[i_prop] = (TYPE) expr_temp.subs_numerical(var);
        }
        
        if (matprop[i].type=="isotrope" and matprop[i].comp=="elastique") {
            matprop[i].coef.push_back(mat_prop_temp[0]);   // E
            matprop[i].coef.push_back(mat_prop_temp[1]);   // nu
            matprop[i].coef.push_back(mat_prop_temp[2]);   // alpha
            matprop[i].coef.push_back(0);                                              // deltaT
        } else if (matprop[i].type=="isotrope" and matprop[i].comp=="visqueux") {
            matprop[i].coef.push_back(mat_prop_temp[0]);   // E
            matprop[i].coef.push_back(mat_prop_temp[1]);   // nu
            matprop[i].coef.push_back(mat_prop_temp[4]);   // viscosite
            matprop[i].coef.push_back(mat_prop_temp[2]);   // alpha
            matprop[i].coef.push_back(0);                                              // deltaT
        } else if (matprop[i].type=="orthotrope") {
            matprop[i].coef.push_back(mat_prop_temp[14]);   // E1
            matprop[i].coef.push_back(mat_prop_temp[15]);   // E2
            matprop[i].coef.push_back(mat_prop_temp[16]);   // E3
            
            matprop[i].coef.push_back(mat_prop_temp[20]);   // nu12
            matprop[i].coef.push_back(mat_prop_temp[22]);   // nu13
            matprop[i].coef.push_back(mat_prop_temp[21]);   // nu23
            
            matprop[i].coef.push_back(mat_prop_temp[17]);   // G12
            matprop[i].coef.push_back(mat_prop_temp[19]);   // G13
            matprop[i].coef.push_back(mat_prop_temp[18]);   // G23

            for(int d=0; d<data_user.dim; d++){
                matprop[i].direction[0][d]=mat_prop_temp[d+5];
                matprop[i].direction[0][d]=mat_prop_temp[d+8];
            }
            
            std::cout << "assignation_materials_sst.h " << matprop[i].direction[0] << " v1 et v2 " << matprop[i].direction[1] << std::endl;
            //normalisation des directions d'orthotropie
            matprop[i].direction[0]=matprop[i].direction[0]/norm_2(matprop[i].direction[0]);
            matprop[i].direction[1]=matprop[i].direction[1]/norm_2(matprop[i].direction[1]);

            //coefficients thermiques
            matprop[i].coefth.resize(4);
            matprop[i].coefth[0]=mat_prop_temp[23];         //alpha_1
            matprop[i].coefth[0]=mat_prop_temp[24];         //alpha_2
            matprop[i].coefth[0]=mat_prop_temp[25];         //alpha_3
            matprop[i].coefth[3]=0;
            
            if (matprop[i].comp=="endommageable") {
                //parametres d'endommagement
                matprop[i].param_damage.Yo = mat_prop_temp[26];
                matprop[i].param_damage.Yop = mat_prop_temp[27];
                matprop[i].param_damage.Ysp = mat_prop_temp[28];
                matprop[i].param_damage.Yc = mat_prop_temp[29];
                matprop[i].param_damage.Ycp = mat_prop_temp[30];
                matprop[i].param_damage.b = mat_prop_temp[31];
            }
        }
    }
};




template<class TV3>
void read_material_properties(TV3 &matprop, Param &process, const XmlNode &n) {

    XmlNode nmat=n.get_element("materials");
    unsigned nbmat = nmat.nb_elements("coefficients");

    matprop.resize(nbmat);
    for(unsigned i=0;i<nbmat;++i) {
        XmlNode nc = nmat.get_element("coefficients",i);
        int id;
        nc.get_attribute( "identificateur", id );
        nc.get_attribute( "type", matprop[id].type );

        string resolution;
        if (TV3::template SubType<0>::T::dim==2) {
            nc.get_attribute( "resolution", resolution );
            if (resolution=="contrainte_plane")
                matprop[id].resolution=1;
            else if (resolution=="deformation_plane")
                matprop[id].resolution=0;
            else {
                std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
                assert(0);
            }
        }
        else
            matprop[id].resolution=0;


        std::vector<Ex> symbols;
        if (TV3::template SubType<0>::T::dim==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (TV3::template SubType<0>::T::dim==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }
        string fvolstr;
        nc.get_attribute( "f_vol", fvolstr );
        Vec<string> vstr=tokenize(fvolstr,';');
        Vec<Ex> expr;
        expr.resize(TV3::template SubType<0>::T::dim);
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(vstr[d2],symbols);
        }
        Vec<double,TV3::template SubType<0>::T::dim> data;
        Ex::MapExNum var;
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= 0.;
        }
        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            data[d2] = (double)expr[d2].subs_numerical(var);

        matprop[id].f_vol=data;
//         std::cout << "Pour le materiau  " << id << " : " << data << std::endl;

        if (matprop[id].type=="isotrope") {
            double E;
            nc.get_attribute( "elastic_modulus", E );
            double nu;
            nc.get_attribute( "poisson_ratio", nu );
            matprop[id].coef.push_back(E);
            matprop[id].coef.push_back(nu);
            double alpha,deltaT;
            nc.get_attribute( "alpha", alpha );
            deltaT=process.properties->deltaT;
            Vec<double> coefth(alpha,deltaT);
            matprop[id].coefth.append(coefth);
            if(process.temps->type_de_calcul=="Qstat")
                matprop[id].dt = process.temps->dt;

        } else if (matprop[id].type=="viscoelas") {
            double E;
            nc.get_attribute( "elastic_modulus", E );
            double nu;
            nc.get_attribute( "poisson_ratio", nu );
            double viscosite;
            nc.get_attribute( "viscosite", viscosite );
            matprop[id].coef.push_back(E);
            matprop[id].coef.push_back(nu);
            matprop[id].coef.push_back(viscosite);
            double alpha,deltaT;
            nc.get_attribute( "alpha", alpha );
            deltaT=process.properties->deltaT;
            Vec<double> coefth(alpha,deltaT);
            matprop[id].coefth.append(coefth);
            if(process.temps->type_de_calcul=="Qstat")
                matprop[id].dt = process.temps->dt;
        } else if (matprop[id].type=="orthotrope" or matprop[id].type=="orthotrope_DPG" or matprop[id].type=="orthotrope_damage") {
            double E1,E2,E3,nu12,nu13,nu23,G12,G13,G23;
            nc.get_attribute( "E1", E1 );
            nc.get_attribute( "E2", E2 );
            nc.get_attribute( "E3", E3 );
            nc.get_attribute( "nu12", nu12 );
            nc.get_attribute( "nu13", nu13 );
            nc.get_attribute( "nu23", nu23 );
            nc.get_attribute( "G12", G12 );
            nc.get_attribute( "G13", G13 );
            nc.get_attribute( "G23", G23);
            Vec<double> coef(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23);
            matprop[id].coef.append(coef);

            string exprv1 = nc.get_attribute( "v1" ),exprv2 = nc.get_attribute( "v2" );
            Vec<string> v1=tokenize(exprv1,';');
            Vec<string> v2=tokenize(exprv2,';');

            Ex Pi = symbol("Pi");
            std::vector<Ex> symbols;
            symbols.push_back(Pi);
            for(unsigned q=0;q<v1.size();++q) {
                Ex res1 = read_ex(v1[q].c_str(),symbols);
                Ex res2 = read_ex(v2[q].c_str(),symbols);

                matprop[id].direction[0][q]=res1.subs_numerical(Pi,M_PI);
                matprop[id].direction[1][q]=res2.subs_numerical(Pi,M_PI);
            }
            std::cout << "assignation_materials_sst.h " << matprop[id].direction[0] << " v1 et v2 " << matprop[id].direction[1] << std::endl;
            //normalisation des directions d'orthotropie
            matprop[id].direction[0]=matprop[id].direction[0]/norm_2(matprop[id].direction[0]);
            matprop[id].direction[1]=matprop[id].direction[1]/norm_2(matprop[id].direction[1]);

            //coefficients thermiques
            matprop[id].coefth.resize(4);
            nc.get_attribute( "alpha_1", matprop[id].coefth[0] );
            nc.get_attribute( "alpha_2", matprop[id].coefth[1] );
            nc.get_attribute( "alpha_3", matprop[id].coefth[2] );

            double deltaT=process.properties->deltaT;
            matprop[id].coefth[3]=deltaT;
            if(process.temps->type_de_calcul=="Qstat")
                matprop[id].dt = process.temps->dt;
            //if (matprop[id].type=="orthotrope_DPG"){
            //   matprop[id].epshorsplan = process.properties->epshorsplan;
            //   }
            if (matprop[id].type=="orthotrope_damage") {
                //parametres d'endommagement
                nc.get_attribute( "Yo", matprop[id].param_damage.Yo );
                nc.get_attribute( "Yop", matprop[id].param_damage.Yop );
                nc.get_attribute( "Ysp", matprop[id].param_damage.Ysp );
                nc.get_attribute( "Yc", matprop[id].param_damage.Yc );
                nc.get_attribute( "Ycp", matprop[id].param_damage.Ycp );
                nc.get_attribute( "b", matprop[id].param_damage.b );
            }

        }
    }
};

/** \ingroup Materiaux
\brief Lecture des proprietes materiau d'interface 
 
On recherche le bloc <proprietes_interfaces>  </ proprietes_interfaces> puis on détermine le nombre de propriétés différentes en considérant le nombre de champs <coefficients />. On crée ainsi un vecteur de propriétés d'interface.
 
On a la possibilité d'entrer différentes propriétés particulières : 
- interface de type contact renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_sst" coeffrottement="0.3" num_sst="0 1" name="cube_cube"/>
\endcode
- interface de type contact renseignée par une boite déterminée par ses deux points extrèmes (classés par ordre croissant).
\code 
<coefficients type="contact_box" coeffrottement="0.3" box="1 0 0 1 10 10" name="cube_cube"/>
\endcode
Le champ box contient donc le point inférieur gauche puis le point supérieur droit (4 composantes en 2d 6 en 3d). Toutes les interfaces dont la boite est incluse dans cette boite sont repérées. Attention pour les interfaces courbes (la boite n'est plus un plan ou une droite).
- interface de contact avec jeu renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_jeu_sst" coeffrottement="0.3" num_sst="0 1" jeu="x*x+y+1" name="cube_cube"/>
\endcode
Le jeu est une fonction quelconque des coordonnées des interfaces concernées dans le repère x, y, z. S'il est rentré comme précédement, il sera considéré comme étant normal aux surfaces. On peut également le rentrer complet <tt>jeu="x*x;0;0"</tt> ce qui construira un champs par point fonction de l'espace.
- interface de contact avec jeu renseignée par une boite déterminée par ses deux points extrêmes (classés par ordre croissant).
\code 
<coefficients type="contact_jeu_box" coeffrottement="0.3" box="0 0 0 0 0 0" jeu="x*x+y+1" name="fibre-matrice"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de type cohésive renseignée par une boite.
\code 
<coefficients type="cohesive" coeffrottement="0.3" kn="0.12" kt="0.12" knc="0.12" gamma="0" alpha="0" Yc="0" Yo="0" n="2.5" name="fibre-matrice"/>
\endcode
- interface de type discrète
\code
<coefficients type="discrete" coeffrottement="0.3" Gcrit="0.12"  name="fibre-matrice"/>
\endcode
Lorsque le taux de restitution critique est inférieur à une valeur critique cette interface est parfaite sinon elle est de type contact avec frottement. Le nom est ici important puisque c'est à partir de celui-ci que les propriétés matériaux sont assignées (on recherche de quel type sont les sous-structures adjacentes)
- interface de type jeu imposé renseignée par le numéro des deux sous-structures adjacentes :
\code 
<coefficients type="jeu_impose_sst" num_sst="0 1" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de contact avec jeu renseignée par une boite déterminée par ses deux points extrêmes (classés par ordre croissant).
\code 
<coefficients type="jeu_impose_box" box="0 0 0 0 0 0" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
Même remarque sur le jeu que pour l'interface contact_jeu_sst.
 
*/

template<class TV4>
void read_propinter(TV4 &propinter,const DataUser &data_user) {
    unsigned nbliaisons = data_user.behaviour_links.size();
    propinter.resize(nbliaisons);
    for(unsigned i=0;i<nbliaisons;++i) {
        propinter[i].id = data_user.behaviour_links[i].id;
        if(data_user.behaviour_links[i].type == "contact"){
            if(data_user.behaviour_links[i].comp_complexe == ""){
                propinter[i].type = "contact_jeu_sst";
                propinter[i].comp="Contact_jeu";
                break;
            }else if(data_user.behaviour_links[i].comp_complexe == "Ca"){
                propinter[i].type = "cohesive";
                propinter[i].comp="Cohesive";
                break;
            }
            break;
        }else if(data_user.behaviour_links[i].type == "parfait"){
            if(data_user.behaviour_links[i].comp_complexe == ""){
                propinter[i].type = "parfait";
                propinter[i].comp="Parfait";
                break;
            }else if(data_user.behaviour_links[i].comp_complexe == "Ca"){
                propinter[i].type = "parfait";
                propinter[i].comp="Parfait";
                break;
            }
            break;
        }
        
        
        std::vector<Ex> symbols;
        if (TV4::template SubType<0>::T::dim==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (TV4::template SubType<0>::T::dim==3) {
            symbols.push_back("x");
            symbols.push_back("y");
            symbols.push_back("z");
        }
        
        Vec< TYPE > link_prop_temp;
        link_prop_temp.resize(data_user.behaviour_links[i].link_prop.size());
        for(int i_prop=0; i_prop<data_user.behaviour_links[i].link_prop.size(); i_prop++){
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_links[i].link_prop[i_prop],symbols);
            Ex::MapExNum var;
            for(unsigned d2=0;d2<TV4::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var[symbols[d2]]= 0.;
            }
            link_prop_temp[i_prop] = (TYPE) expr_temp.subs_numerical(var);
        }
        
        propinter[i].coeffrottement = link_prop_temp[0];                                // coeff frottement
        propinter[i].jeu = data_user.behaviour_links[i].link_prop[1];                   // jeux ou epaisseur negative
        propinter[i].Gcrit = link_prop_temp[7];                                         // limite en rupture    
    }
}


template<class TV4>
void read_propinter(TV4 &propinter,const XmlNode &n) {

    XmlNode nmat=n.get_element("proprietes_interfaces");
    unsigned nbmat = nmat.nb_elements("coefficients");

    propinter.resize(nbmat);

    for(unsigned i=0;i<nbmat;++i) {

        XmlNode nc = nmat.get_element("coefficients",i);
        nc.get_attribute( "type", propinter[i].type);
        nc.get_attribute( "name",propinter[i].name);
        if (propinter[i].type=="contact_sst") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            nc.get_attribute( "num_sst", propinter[i].num_sst);
            propinter[i].box[0].set(0.);
            propinter[i].box[1].set(0.);
            propinter[i].comp="Contact";
        } else if (propinter[i].type=="contact_box") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            Vec<typename TV4::template SubType<0>::T::T> box;
            nc.get_attribute( "box",box);
            propinter[i].box[0]=box[range((int)TV4::template SubType<0>::T::dim)]
                                ;
            propinter[i].box[1]=box[range((int)TV4::template SubType<0>::T::dim,(int)(2*TV4::template SubType<0>::T::dim))]
                                ;
            propinter[i].num_sst.set(0);
            propinter[i].comp="Contact";
        } else if (propinter[i].type=="contact_jeu_sst") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            nc.get_attribute( "num_sst", propinter[i].num_sst);
            propinter[i].box[0].set(0.);
            propinter[i].box[1].set(0.);
            nc.get_attribute( "jeu",propinter[i].jeu);
            propinter[i].comp="Contact_jeu";
        } else if (propinter[i].type=="contact_jeu_box") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            Vec<typename TV4::template SubType<0>::T::T> box;
            nc.get_attribute( "box",box);
            propinter[i].box[0]=box[range((int)TV4::template SubType<0>::T::dim)]
                                ;
            propinter[i].box[1]=box[range((int)TV4::template SubType<0>::T::dim,(int)(2*TV4::template SubType<0>::T::dim))]
                                ;
            propinter[i].num_sst.set(0);
            nc.get_attribute( "jeu",propinter[i].jeu);
            propinter[i].comp="Contact_jeu";
        } else if (propinter[i].type=="contact_jeu_physique") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            nc.get_attribute( "num_sst", propinter[i].num_sst);
            propinter[i].box[0].set(0.);
            propinter[i].box[1].set(0.);
            propinter[i].comp="Contact_jeu_physique";
        } else if (propinter[i].type=="discrete") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            nc.get_attribute( "Gcrit", propinter[i].Gcrit);
            Vec<typename TV4::template SubType<0>::T::T> box;
            nc.get_attribute( "box",box);
            propinter[i].box[0]=box[range((int)TV4::template SubType<0>::T::dim)]
                                ;
            propinter[i].box[1]=box[range((int)TV4::template SubType<0>::T::dim,(int)(2*TV4::template SubType<0>::T::dim))]
                                ;
            propinter[i].num_sst.set(0);
            propinter[i].comp="Parfait";
        } else if (propinter[i].type=="cohesive") {
            nc.get_attribute( "coeffrottement", propinter[i].coeffrottement);
            Vec<typename TV4::template SubType<0>::T::T> box;
            nc.get_attribute( "box",box);
            propinter[i].box[0]=box[range((int)TV4::template SubType<0>::T::dim)]
                                ;
            propinter[i].box[1]=box[range((int)TV4::template SubType<0>::T::dim,(int)(2*TV4::template SubType<0>::T::dim))]
                                ;
            propinter[i].num_sst.set(0);
            nc.get_attribute( "kn",propinter[i].param_damage.kn);
            nc.get_attribute( "knc",propinter[i].param_damage.knc);
            nc.get_attribute( "kt",propinter[i].param_damage.kt);
            nc.get_attribute( "alpha",propinter[i].param_damage.alpha);
            nc.get_attribute( "gamma",propinter[i].param_damage.gamma);
            nc.get_attribute( "Yc",propinter[i].param_damage.Yc);
            nc.get_attribute( "Yo",propinter[i].param_damage.Yo);
            nc.get_attribute( "n",propinter[i].param_damage.n);
            propinter[i].comp="Cohesive";
        } else if (propinter[i].type=="jeu_impose_sst") {
            propinter[i].coeffrottement=0;
            nc.get_attribute( "num_sst", propinter[i].num_sst);
            propinter[i].box[0].set(0.);
            propinter[i].box[1].set(0.);
            nc.get_attribute( "jeu",propinter[i].jeu);
            nc.get_attribute( "nbpas",propinter[i].nbpastempsimpos,1);
            propinter[i].comp="Jeu_impose";
        } else if (propinter[i].type=="jeu_impose_box") {
            propinter[i].coeffrottement=0;
            Vec<typename TV4::template SubType<0>::T::T> box;
            nc.get_attribute( "box",box);
            propinter[i].box[0]=box[range((int)TV4::template SubType<0>::T::dim)]
                                ;
            propinter[i].box[1]=box[range((int)TV4::template SubType<0>::T::dim,(int)(2*TV4::template SubType<0>::T::dim))]
                                ;
            propinter[i].num_sst.set(0);
            nc.get_attribute( "jeu",propinter[i].jeu);
            nc.get_attribute( "nbpas",propinter[i].nbpastempsimpos,1);
            propinter[i].comp="Jeu_impose";
        } else {
            std::cout << "comportement d'interfaces non implemente : choix : contact_sst, contact_box, contact_jeu_sst, contact_jeu_box, contact_jeu_physique, discrete, cohesive, jeu_impose_sst, jeu_impose_box" << std::endl;
            assert(0);
        }

    }

};

