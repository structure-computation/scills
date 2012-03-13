
#include "codegen/codegen.h"
#include "containers/basicops.h"
#include "DataUser.h"
#include "../../../SC_code/src/COMPUTE/FieldStructureUser.h"

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
void read_material_properties(TV3 &matprop, Param &process, DataUser &data_user , BasicVec<BasicVec<TYPE> > &mat_prop_temp) {

    unsigned nbmat = data_user.behaviour_materials.size();
//     PRINT(nbmat);
    matprop.resize(nbmat);
    mat_prop_temp.resize(nbmat);
    for(unsigned i=0;i<nbmat;++i) {
        matprop[i].id = data_user.behaviour_materials[i].id;
        matprop[i].type_num = data_user.behaviour_materials[i].type_num;
        matprop[i].type = data_user.behaviour_materials[i].type;
        matprop[i].comp = data_user.behaviour_materials[i].comp;
        if(data_user.dim == 2){
            if (data_user.behaviour_materials[i].resolution =="CP")
                matprop[i].resolution=1;
            else if (data_user.behaviour_materials[i].resolution =="DP")
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
        
        if(data_user.options.Multiresolution_on==1){
            //ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                String var="V"; var<<i_par; //nom de la variable de multiresolution
                symbols.push_back(var.c_str());
            }
        }

        Ex Pi = symbol("PI");
        symbols.push_back(Pi);
        
        Vec<std::string> vstr;
        vstr.resize(data_user.dim, "0");
        
// 	std::cout << "data_user.behaviour_bc_volume[1].select = " << data_user.behaviour_bc_volume[1].select << std::endl;
// 	std::cout << "data_user.behaviour_bc_volume[0].select = " << data_user.behaviour_bc_volume[0].select << std::endl;
        for(int i_fvol=0; i_fvol<data_user.behaviour_bc_volume.size(); i_fvol++){
            if(data_user.behaviour_bc_volume[i_fvol].select){
                for(int d=0; d<data_user.dim; d++){
                    vstr[d] += " + " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[d] + " * " + data_user.behaviour_bc_volume[i_fvol].step[0].CLv_step_prop[6];
                }
            }
        }
        std::cout << "force volumique 0 = " << vstr[0] << std::endl;
        std::cout << "force volumique 1 = " << vstr[1] << std::endl;
        std::cout << "force volumique 2 = " << vstr[2] << std::endl;
        
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
        if(data_user.options.Multiresolution_on==1){
            //evaluation des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var[symbols[TV3::template SubType<0>::T::dim+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        var[symbols[TV3::template SubType<0>::T::dim+data_user.Multiresolution_parameters.size()]]=M_PI;

        for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            data[d2] = (double)expr[d2].subs_numerical(var);

        matprop[i].f_vol_e=vstr;
        matprop[i].f_vol=data;
	
        
        mat_prop_temp[i].resize(data_user.behaviour_materials[i].mat_prop.size());
//         PRINT(mat_prop_temp.size());
        for(int i_prop=0; i_prop<data_user.behaviour_materials[i].mat_prop.size(); i_prop++){
            if(data_user.behaviour_materials[i].mat_prop[i_prop] == ""){
                data_user.behaviour_materials[i].mat_prop[i_prop] = "0";
            }
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_materials[i].mat_prop[i_prop],symbols);
            Ex::MapExNum var_temp;
            for(unsigned d2=0;d2<TV3::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var_temp[symbols[d2]]= 0.;
            }
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution aux symboles
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[TV3::template SubType<0>::T::dim+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            var_temp[symbols[TV3::template SubType<0>::T::dim+data_user.Multiresolution_parameters.size()]]=M_PI;

            mat_prop_temp[i][i_prop] = (TYPE) expr_temp.subs_numerical(var_temp);
/*            std::cout << "Pour la propriete  " << i_prop << " : " << data << std::endl;*/
        }
/*        std::cout << "Pour le materiau  " << data_user.behaviour_materials[i].id << " : " << data << std::endl;*/
        matprop[i].density = mat_prop_temp[i][3];
        
        std::cout << (matprop[i].type) <<endl;
        if(matprop[i].type.find("isotrope")<matprop[i].type.size()) {                 // comportement isotrope elastique
//            PRINT("comportement isotrope elastique");
            matprop[i].coef.push_back(mat_prop_temp[i][0]);   /// E
            matprop[i].coef.push_back(mat_prop_temp[i][1]);   /// nu
            matprop[i].coefth.push_back(mat_prop_temp[i][2]);   /// alpha
            matprop[i].coefth.push_back(0);                                              /// deltaT
        } else if (matprop[i].comp.find("visqueux")<matprop[i].comp.size() and matprop[i].type.find("isotrope")<matprop[i].type.size()) {          // comportement isotrope elastique visqueux
//             PRINT("comportement isotrope visqueux");
            matprop[i].coef.push_back(mat_prop_temp[i][0]);   /// E
            matprop[i].coef.push_back(mat_prop_temp[i][1]);   /// nu
            matprop[i].coef.push_back(mat_prop_temp[i][4]);   /// viscosite
            matprop[i].coefth.push_back(mat_prop_temp[i][2]);   /// alpha
            matprop[i].coefth.push_back(0);                                              /// deltaT
        } else if (matprop[i].type.find("orthotrope")<matprop[i].type.size()) {          /// orthotrope
            PRINT("comportement orthotrope");
            matprop[i].coef.push_back(mat_prop_temp[i][14]);   /// E1
            matprop[i].coef.push_back(mat_prop_temp[i][15]);   /// E2
            matprop[i].coef.push_back(mat_prop_temp[i][16]);   /// E3
            
            matprop[i].coef.push_back(mat_prop_temp[i][20]);   /// nu12
            matprop[i].coef.push_back(mat_prop_temp[i][22]);   /// nu13
            matprop[i].coef.push_back(mat_prop_temp[i][21]);   /// nu23
            
            matprop[i].coef.push_back(mat_prop_temp[i][17]);   /// G12
            matprop[i].coef.push_back(mat_prop_temp[i][19]);   /// G13
            matprop[i].coef.push_back(mat_prop_temp[i][18]);   /// G23

            for(int d=0; d<data_user.dim; d++){
                matprop[i].direction[0][d]=mat_prop_temp[i][d+5];
                matprop[i].direction[1][d]=mat_prop_temp[i][d+8];
            }
            
            std::cout << "assignation_materials_sst.h " << matprop[i].direction[0] << " v1 et v2 " << matprop[i].direction[1] << std::endl;
            //normalisation des directions d'orthotropie
            matprop[i].direction[0]=matprop[i].direction[0]/norm_2(matprop[i].direction[0]);
            matprop[i].direction[1]=matprop[i].direction[1]/norm_2(matprop[i].direction[1]);

            //coefficients thermiques
            matprop[i].coefth.resize(4);
            matprop[i].coefth[0]=mat_prop_temp[i][23];      ///alpha_1
            matprop[i].coefth[1]=mat_prop_temp[i][24];      ///alpha_2
            matprop[i].coefth[2]=mat_prop_temp[i][25];      ///alpha_3
            matprop[i].coefth[3]=0;
            
            PRINT("comportement plastique");
            //parametres de plasticite
            matprop[i].coefp[0] = mat_prop_temp[i][26];     /// k_p
            matprop[i].coefp[1] = mat_prop_temp[i][27];     /// m_p
            matprop[i].coefp[2] = mat_prop_temp[i][28];     /// R0
            matprop[i].coefp[3] = mat_prop_temp[i][29];     /// couplage
            
            PRINT("comportement endommageable");
            //parametres d'endommagement
            matprop[i].coefendom[0] = mat_prop_temp[i][30];     /// Yo
            matprop[i].coefendom[1] = mat_prop_temp[i][31];     /// Yc
            matprop[i].coefendom[2] = mat_prop_temp[i][32];     /// Ycf
            matprop[i].coefendom[3] = mat_prop_temp[i][33];     /// dmax
            matprop[i].coefendom[4] = mat_prop_temp[i][34];     /// b_c
            matprop[i].coefendom[5] = mat_prop_temp[i][35];     /// effet_retard
            matprop[i].coefendom[6] = mat_prop_temp[i][36];     /// a
            matprop[i].coefendom[7] = mat_prop_temp[i][37];     /// tau_c
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
void read_propinter(TV4 &propinter,const DataUser &data_user, BasicVec<BasicVec<TYPE> > &link_prop_temp) {
    unsigned nbliaisons = data_user.behaviour_links.size();
    propinter.resize(nbliaisons);
    link_prop_temp.resize(nbliaisons);
    for(unsigned i=0;i<nbliaisons;++i) {
        propinter[i].id = data_user.behaviour_links[i].id;
//         PRINT(data_user.behaviour_links[i].type_num);
        if(data_user.behaviour_links[i].type_num == 0){   //parfaite
            propinter[i].type = "parfait";
            propinter[i].comp="Parfait";
        }else if(data_user.behaviour_links[i].type_num == 1){   //elastique
            propinter[i].type = "elastique";
            propinter[i].comp="Parfait";
        }else if(data_user.behaviour_links[i].type_num == 2){
//             propinter[i].type = "contact_jeu_sst";
//             propinter[i].comp="Contact_jeu"; 
            propinter[i].type = "contact_ep";
            propinter[i].comp="Contact_ep";
        }else if(data_user.behaviour_links[i].type_num == 3){
            propinter[i].type = "cohesive";
            propinter[i].comp="Cohesive";
        }else if(data_user.behaviour_links[i].type_num == 4){
            propinter[i].type = "breakable";
            propinter[i].comp="Breakable";
        }else{
            std::cout  << "comportement d'interface non reconnu" << std::endl;  assert(0);
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
        
        if(data_user.options.Multiresolution_on==1){
            //ajout des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++){
                String var="V"; var<<i_par; //nom de la variable de multiresolution
                symbols.push_back(var.c_str());
            }
        }
        
        link_prop_temp[i].resize(data_user.behaviour_links[i].link_prop.size());
        for(int i_prop=0; i_prop<data_user.behaviour_links[i].link_prop.size(); i_prop++){
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_links[i].link_prop[i_prop],symbols);
            Ex::MapExNum var_temp;
            for(unsigned d2=0;d2<TV4::template SubType<0>::T::dim;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var_temp[symbols[d2]]= 0.;
            }
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[TV4::template SubType<0>::T::dim+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            link_prop_temp[i][i_prop] = (TYPE) expr_temp.subs_numerical(var_temp);           
        }
        propinter[i].f_coeffrottement = data_user.behaviour_links[i].link_prop[0];      /// coeff frottement analytique
        propinter[i].coeffrottement = link_prop_temp[i][0];                             /// coeff frottement
        propinter[i].jeu = data_user.behaviour_links[i].link_prop[1];                   /// jeux ou epaisseur negative        
        propinter[i].Gcrit = link_prop_temp[i][7];                                      /// limite en rupture    
        propinter[i].f_R = data_user.behaviour_links[i].link_prop[3];                   /// coeff frottement analytique
        PRINT(link_prop_temp[i]);
    }        

}
