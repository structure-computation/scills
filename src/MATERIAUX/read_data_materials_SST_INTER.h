
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
 
On d�termine ensuite le nombre de blocs <coefficients /> pour cr�er le vecteur de propri�t�s mat�riau (classe SstCarac).
 
L'"identificateur" permet de positionner le mat�riau lu dans le vecteur de propri�t�s. Le num�ro affect� aux sous-structures correspond alors � cet identificateur.
 
Selon la dimension le champ "resolution" est lu ou non. Celui ci peut prendre les valeurs
\code resolution="deformation_plane" ou resolution="contrainte_plane" \endcode
 
Ensuite le champ "type" renseigne sur le comportement de la sous-structure. Selon la valeur obtenue, les diff�rentes propri�t�s sont extraites du fichier xml (les champs n�cessaires sont consultables dans les fichiers de formulations correspondants).
- \code type="isotrope" \endcode
- \code type="orthotrope" \endcode
- \code type="viscoelas" \endcode
- \code type="orthotrope_damage" \endcode
 
On aura donc les possibilit�s suivantes : 
-   On peut choisir un mat�riau isotrope  :
   \code 
   <coefficients type="isotrope" resolution="contrainte_plane" identificateur="0" name="fibre"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6"/>
   \endcode
-   ou bien un mat�riau orthotrope
   \code
   <coefficients type="orthotrope" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0"/>
   \endcode
-   ou encore un mat�riau orthotrope endommageable pour lequel on ajoute 
   \code
   <coefficients type="orthotrope_damage" resolution="contrainte_plane" identificateur="0" E1="157e3" E2="8500" E3="8500" nu12="0.29" nu13="0.29" nu23="0.4" G12="5000" G13="5000" G23="3000" unit="MPa" thickness="1" alpha_1="2.3e-6" alpha_2="30e-6" alpha_3="30e-6" v1="0;1;0" v2="-1;0;0" Yo="0.0961" Yop="0.0961" Ycp="10.8" Yc="10.8" b="2.5" Ysp="0.70"/>
   \endcode
-   ou encore un mat�riau viscoelastique (orthotrope ou isotrope)
   \code
   <coefficients type="isotrope" resolution="contrainte_plane" identificateur="0"  elastic_modulus="200e3" poisson_ratio="0.3" unit="MPa" thickness="1" alpha="2e-6" viscosite="0.1" />
   \endcode
*/




template<class TV3>
void read_material_properties(TV3 &matprops, Param &process, DataUser &data_user , BasicVec<BasicVec<TYPEREEL> > &mat_prop_temp) {

    unsigned nbmat = data_user.behaviour_materials.size();
//     PRINT(nbmat);
    matprops.resize(nbmat);
    mat_prop_temp.resize(nbmat);
    for(unsigned i=0;i<nbmat;++i) {
        matprops[i].id = data_user.behaviour_materials[i].id;
        matprops[i].type_num = data_user.behaviour_materials[i].type_num;
        matprops[i].type = data_user.behaviour_materials[i].type;
        matprops[i].comp = data_user.behaviour_materials[i].comp;
        if(data_user.dim == 2){
            if (data_user.behaviour_materials[i].resolution =="CP")
                matprops[i].resolution=1;
            else if (data_user.behaviour_materials[i].resolution =="DP")
                matprops[i].resolution=0;
            else {
                std::cout << "type de resolution non implemente : choix contrainte_plane ou deformation_plane" << std::endl;
                assert(0);
            }
        }else{
            matprops[i].resolution=0;
        }

        std::vector<Ex> symbols;
        if (DIM==2) {
            symbols.push_back("x");
            symbols.push_back("y");
        }
        else if (DIM==3) {
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
        expr.resize(DIM);
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            expr[d2] = read_ex(vstr[d2],symbols);
        }
        
        Vec<double,DIM> data;
        Ex::MapExNum var;
        for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
            var[symbols[d2]]= 0.;
        }
        if(data_user.options.Multiresolution_on==1){
            //evaluation des variables de multiresolution aux symboles
            for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                var[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
        }
        var[symbols[DIM+data_user.Multiresolution_parameters.size()]]=M_PI;

        for(unsigned d2=0;d2<DIM;++d2)//boucle sur les inconnues possibles (dimension des vecteurs)
            data[d2] = (double)expr[d2].subs_numerical(var);

        matprops[i].f_vol_e=vstr;
        matprops[i].f_vol=data;
	
        
        mat_prop_temp[i].resize(data_user.behaviour_materials[i].mat_prop.size());
//         PRINT(mat_prop_temp.size());
        for(int i_prop=0; i_prop<data_user.behaviour_materials[i].mat_prop.size(); i_prop++){
            if(data_user.behaviour_materials[i].mat_prop[i_prop] == ""){
                data_user.behaviour_materials[i].mat_prop[i_prop] = "0";
            }
            Ex expr_temp;
            expr_temp = read_ex(data_user.behaviour_materials[i].mat_prop[i_prop],symbols);
            Ex::MapExNum var_temp;
            for(unsigned d2=0;d2<DIM;++d2) {//boucle sur les inconnues possibles (dimension des vecteurs)
                var_temp[symbols[d2]]= 0.;
            }
            if(data_user.options.Multiresolution_on==1){
                //evaluation des variables de multiresolution aux symboles
                for(unsigned i_par=0;i_par< data_user.Multiresolution_parameters.size() ;i_par++)
                    var_temp[symbols[DIM+i_par]]=data_user.Multiresolution_parameters[i_par].current_value;
            }
            var_temp[symbols[DIM+data_user.Multiresolution_parameters.size()]]=M_PI;

            mat_prop_temp[i][i_prop] = (TYPE) expr_temp.subs_numerical(var_temp);
/*            std::cout << "Pour la propriete  " << i_prop << " : " << data << std::endl;*/
        }
/*        std::cout << "Pour le materiau  " << data_user.behaviour_materials[i].id << " : " << data << std::endl;*/
        matprops[i].density = mat_prop_temp[i][3];
        
        std::cout << (matprops[i].type) <<endl;
        if(matprops[i].type.find("isotrope")<matprops[i].type.size()) {/// comportement isotrope elastique
//            PRINT("comportement isotrope elastique");
            matprops[i].elastic_modulus = mat_prop_temp[i][0];   /// E
            matprops[i].poisson_ratio   = mat_prop_temp[i][1];   /// nu
            matprops[i].alpha           = mat_prop_temp[i][2];   /// alpha
            matprops[i].deltaT          = 0;                     /// deltaT
        } else if (matprops[i].comp.find("visqueux")<matprops[i].comp.size() and matprops[i].type.find("isotrope")<matprops[i].type.size()) {/// comportement isotrope elastique visqueux
//             PRINT("comportement isotrope visqueux");
            matprops[i].elastic_modulus = mat_prop_temp[i][0];   /// E
            matprops[i].poisson_ratio   = mat_prop_temp[i][1];   /// nu
            matprops[i].alpha           = mat_prop_temp[i][2];   /// alpha
            matprops[i].deltaT          = 0;                     /// deltaT
            matprops[i].viscosite       = mat_prop_temp[i][4];   /// viscosite
        } else if (matprops[i].type.find("orthotrope")<matprops[i].type.size() or
                   matprops[i].type.find("mesomodele")<matprops[i].type.size()) {/// orthotrope
            PRINT("comportement orthotrope");
            matprops[i].elastic_modulus_1 = mat_prop_temp[i][14];   /// E1
            matprops[i].elastic_modulus_2 = mat_prop_temp[i][15];   /// E2
            matprops[i].elastic_modulus_3 = mat_prop_temp[i][16];   /// E3
            
            matprops[i].poisson_ratio_12 = mat_prop_temp[i][20];   /// nu12
            matprops[i].poisson_ratio_13 = mat_prop_temp[i][22];   /// nu13
            matprops[i].poisson_ratio_23 = mat_prop_temp[i][21];   /// nu23
            
            matprops[i].shear_modulus_12 = mat_prop_temp[i][17];   /// G12
            matprops[i].shear_modulus_13 = mat_prop_temp[i][19];   /// G13
            matprops[i].shear_modulus_23 = mat_prop_temp[i][18];   /// G23

            for(int d=0; d<data_user.dim; d++){
                matprops[i].v1[d]=mat_prop_temp[i][d+5];
                matprops[i].v2[d]=mat_prop_temp[i][d+8];
            }
            
            std::cout << "assignation_materials_sst.h " << matprops[i].v1 << " v1 et v2 " << matprops[i].v2 << std::endl;
            //normalisation des directions d'orthotropie
            matprops[i].v1=matprops[i].v1/norm_2(matprops[i].v1);
            matprops[i].v2=matprops[i].v2/norm_2(matprops[i].v2);

            //coefficients thermiques
            matprops[i].alpha_1 = mat_prop_temp[i][23];      ///alpha_1
            matprops[i].alpha_2 = mat_prop_temp[i][24];      ///alpha_2
            matprops[i].alpha_3 = mat_prop_temp[i][25];      ///alpha_3
            matprops[i].deltaT  = 0;                         /// deltaT
            
            PRINT("comportement plastique");
            //parametres de plasticite
            matprops[i].k_p              = mat_prop_temp[i][26];     /// k_p
            matprops[i].m_p              = mat_prop_temp[i][27];     /// m_p
            matprops[i].R0               = mat_prop_temp[i][28];     /// R0
            matprops[i].coefvm_composite = mat_prop_temp[i][29];     /// couplage
            
            PRINT("comportement endommageable");
            //parametres d'endommagement
            matprops[i].Yo           = mat_prop_temp[i][30];     /// Yo
            matprops[i].Yc           = mat_prop_temp[i][31];     /// Yc
            matprops[i].Ycf          = mat_prop_temp[i][32];     /// Ycf
            matprops[i].dmax         = mat_prop_temp[i][33];     /// dmax
            matprops[i].b_c          = mat_prop_temp[i][34];     /// b_c
            matprops[i].effet_retard = mat_prop_temp[i][35];     /// effet_retard
            matprops[i].a            = mat_prop_temp[i][36];     /// a
            matprops[i].tau_c        = mat_prop_temp[i][37];     /// tau_c
        }
    }
};


/** \ingroup Materiaux
\brief Lecture des proprietes materiau d'interface 
 
On recherche le bloc <proprietes_interfaces>  </ proprietes_interfaces> puis on d�termine le nombre de propri�t�s diff�rentes en consid�rant le nombre de champs <coefficients />. On cr�e ainsi un vecteur de propri�t�s d'interface.
 
On a la possibilit� d'entrer diff�rentes propri�t�s particuli�res : 
- interface de type contact renseign�e par le num�ro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_sst" coeffrottement="0.3" num_sst="0 1" name="cube_cube"/>
\endcode
- interface de type contact renseign�e par une boite d�termin�e par ses deux points extr�mes (class�s par ordre croissant).
\code 
<coefficients type="contact_box" coeffrottement="0.3" box="1 0 0 1 10 10" name="cube_cube"/>
\endcode
Le champ box contient donc le point inf�rieur gauche puis le point sup�rieur droit (4 composantes en 2d 6 en 3d). Toutes les interfaces dont la boite est incluse dans cette boite sont rep�r�es. Attention pour les interfaces courbes (la boite n'est plus un plan ou une droite).
- interface de contact avec jeu renseign�e par le num�ro des deux sous-structures adjacentes :
\code 
<coefficients type="contact_jeu_sst" coeffrottement="0.3" num_sst="0 1" jeu="x*x+y+1" name="cube_cube"/>
\endcode
Le jeu est une fonction quelconque des coordonn�es des interfaces concern�es dans le rep�re x, y, z. S'il est rentr� comme pr�c�dement, il sera consid�r� comme �tant normal aux surfaces. On peut �galement le rentrer complet <tt>jeu="x*x;0;0"</tt> ce qui construira un champs par point fonction de l'espace.
- interface de contact avec jeu renseign�e par une boite d�termin�e par ses deux points extr�mes (class�s par ordre croissant).
\code 
<coefficients type="contact_jeu_box" coeffrottement="0.3" box="0 0 0 0 0 0" jeu="x*x+y+1" name="fibre-matrice"/>
\endcode
M�me remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de type coh�sive renseign�e par une boite.
\code 
<coefficients type="cohesive" coeffrottement="0.3" kn="0.12" kt="0.12" knc="0.12" gamma="0" alpha="0" Yc="0" Yo="0" n="2.5" name="fibre-matrice"/>
\endcode
- interface de type discr�te
\code
<coefficients type="discrete" coeffrottement="0.3" Gcrit="0.12"  name="fibre-matrice"/>
\endcode
Lorsque le taux de restitution critique est inf�rieur � une valeur critique cette interface est parfaite sinon elle est de type contact avec frottement. Le nom est ici important puisque c'est � partir de celui-ci que les propri�t�s mat�riaux sont assign�es (on recherche de quel type sont les sous-structures adjacentes)
- interface de type jeu impos� renseign�e par le num�ro des deux sous-structures adjacentes :
\code 
<coefficients type="jeu_impose_sst" num_sst="0 1" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
M�me remarque sur le jeu que pour l'interface contact_jeu_sst.
- interface de contact avec jeu renseign�e par une boite d�termin�e par ses deux points extr�mes (class�s par ordre croissant).
\code 
<coefficients type="jeu_impose_box" box="0 0 0 0 0 0" jeu="x*x+y+1" name="vis-ecrou"/>
\endcode
M�me remarque sur le jeu que pour l'interface contact_jeu_sst.
 
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
