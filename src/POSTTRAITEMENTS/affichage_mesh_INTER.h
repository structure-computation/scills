#ifndef AFFICH_MESH_INTER_H
#define AFFICH_MESH_INTER_H

using namespace LMT;
using namespace std;

// Eclate des interfaces

struct eclat_INTER_struct {
   template<class INTER> void operator()(INTER &Inter, typename INTER::T ecl) const{
   typedef Vec<typename INTER::T,INTER::dim> TV;
   // calcul cdg du maillage
   
   TV G=barycenter_constant_rho(*Inter.side[0].mesh);
   // translation du maillage
   TV trans=G*ecl;
   // modification du champ qtrans
   apply(Inter.side[0].mesh->node_list,modif_qtrans(),trans);
   }
};


/** \ingroup  Post_Traitement 
\brief Eclaté des maillages d'interface
*/
template<class TV2> void eclat_INTER(TV2 &Inter,typename TV2::template SubType<0>::T::T ecl) {
   apply(Inter,eclat_INTER_struct(),ecl);
};


/** \ingroup  Post_Traitement 
\brief Application du type d'interface et d'un numéro à chaque élément d'interface
*/
struct apply_type_elem_interface{
//template<class TE> void operator() ( TE &e, unsigned i,const int type,const int num, int &ii) const{
template<class TE> void operator() ( TE &e,const int type,const int num, int &ii) const{
    e.type=type;
    e.num=num;
    e.d=ii;
    ii+=1;
  }
};



//*****************************
// affichage des INTERFACES
//********************************
/** \ingroup  Post_Traitement 
\brief Affichage du type de l'interface (CL,...) et du numero de l'interface

Selon le type d'interface un numéro est assigné à chaque élément d'interface : 
- 0 : interface de type déplacement
- 1 : interface de type effort
- 2 : interface de type symetrie 
- 3 : interface de type deplacement normal
- 4 : interface de type parfait
- 5 : interface de type contact ou contact avec jeu ou contact physique
- 6 : interface de type jeu imposé
- 7 : interface périodique
- 8 : interface d'un autre type
On visualise ensuite le type d'interface par le champ "type" dans paraview.

Le numéro de l'interface est aussi assigné au champ "num"
*/
template<class TV2, class TV1> void affich_INTER(TV2 &Inter,TV1 &S, Param &process) {
    string typemail=process.affichage->type_affichage;
    string save=process.affichage->save;
    unsigned data=process.affichage->side;
    string nom_generique = process.affichage->repertoire_save +typemail;

    system(("mkdir -p "+process.affichage->repertoire_save).c_str());

    ostringstream ss;
    ss<<nom_generique << "_"<<process.rank<< "_";
    string namefile(ss.str());
    
    
   //eclate des maillages d'interfaces
   double ecl=1.0;
   eclat_INTER(Inter,ecl);
  
   
   //assignation du numero et du type d'interface
   for(unsigned q=0;q<Inter.size();++q){
                 
         int type=0;
         if (Inter[q].type=="Ext" and Inter[q].comp=="depl"){type=0;}
         else if (Inter[q].type=="Ext" and Inter[q].comp=="effort"){type=1;}
         else if (Inter[q].type=="Ext" and ( Inter[q].comp=="sym" )){type=2;}
         else if (Inter[q].type=="Ext" and ( Inter[q].comp=="depl_normal")){type=3;}
         else if (Inter[q].type=="Int" and Inter[q].comp=="Parfait"){type=4;}
         else if (Inter[q].type=="Int" and (Inter[q].comp=="Contact" or Inter[q].comp=="Contact_jeu" or Inter[q].comp=="Contact_jeu_physique") ){type=5;}
         else if (Inter[q].type=="Int" and Inter[q].comp=="Jeu_impose"){type=6;}
         else if (Inter[q].type=="Ext" and Inter[q].comp=="periodique"){type=7;}
         else {type=8;}
         int num=Inter[q].num;
         int numelem=0;
         apply(Inter[q].side[data].mesh->elem_list,apply_type_elem_interface(),type,num,numelem);
         numelem=0;
         if ( Inter[q].comp=="Contact_jeu_physique" or Inter[q].comp=="periodique") apply(Inter[q].side[1-data].mesh->elem_list,apply_type_elem_interface(),type,num,numelem);
      }

   //affichage
   DisplayParaview dp;
   typename TV2::template SubType<0>::T::TMESH meshglob;
    for(unsigned i=0;i<Inter.size();++i) {
        if (S[Inter[i].vois[data*2]].num_proc==process.rank){
                meshglob.append(*Inter[i].side[data].mesh);
                if ( Inter[i].comp=="Contact_jeu_physique" or Inter[i].comp=="periodique")
                   meshglob.append(*Inter[i].side[1-data].mesh);
        }
    }
   dp.add_mesh(meshglob,namefile,Vec<string>("num","type","qtrans","d"));

   if(save=="display" and process.size==1 and process.affichage->command_file=="")
      dp.exec();
};


// template<class INTER> void affich_INTER_endommagement(Vec<INTER> &Inter, Param &process) {
//    string save = process.affichage->save;
//    for(unsigned q=0;q<Inter.size();++q){
//       int type=0;
//       if (Inter[q].type=="Ext" && Inter[q].comp=="depl"){type=0;}
//       else if (Inter[q].type=="Ext" && Inter[q].comp=="effort"){type=1;}
//       else if (Inter[q].type=="Ext" && Inter[q].comp=="sym"){type=2;}
//       else if (Inter[q].type=="Int" && Inter[q].comp=="Parfait"){type=3;}
//       else if (Inter[q].type=="Int" && Inter[q].comp=="Contact"){type=4;}
//       else if (Inter[q].type=="Int" && Inter[q].comp=="Cohesif"){type=5;}
//       unsigned cpt=0;
//       apply(Inter[q].side[0].mesh.elem_list,apply_type_endommagement_elem_interface(),type,cpt,*Inter[q].parammicro);
//    }
// 
//    DisplayParaview dp;
//    typename INTER::TMESH meshglob;
//    for(unsigned i=0;i<Inter.size();++i)
//       meshglob.append(Inter[i].side[0].mesh);
//            
//    dp.add_mesh(meshglob);
// 
//    if(save=="display")
//       dp.exec();
// 
// };
// 

#endif //AFFICH_MESH_INTER_H
