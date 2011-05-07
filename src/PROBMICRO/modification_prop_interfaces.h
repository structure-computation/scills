


template<class MICRO, class T> void modif_boxdegrad_scale(MICRO &parammicro,T &scale){ 
   parammicro.boxdegrad=scale*parammicro.boxdegrad;
}


template<class TV1, class TV2,class MICRO, class TV4> void modif_inter_degrad(TV2 &Inter,TV1 &S, TV4 &propinter, MICRO &parammicro){
   for(unsigned q=0;q<Inter.size();++q){
      Inter[q].parammicro->degradable=0;
      if (Inter[q].type=="Int"){
            if (inCL(Inter[q].box,parammicro.boxdegrad)==1){
               // differentiation interface fibre-matrice = 0, interface matrice-matrice = 1
               unsigned matvois1=S[Inter[q].vois[0]].typmat,matvois2=S[Inter[q].vois[2]].typmat;
               if (matvois1==matvois2 ) { // matrice-matrice ou //fibre-fibre
                  if (matvois1==0){ //matrice
                  for(unsigned j=0;j<propinter.size();++j){
                     if (propinter[j].name=="matrice-matrice"){
                        Inter[q].parammicro->Gcrit=propinter[j].Gcrit;
                        Inter[q].coefcontact=propinter[j].coeffrottement;
                        Inter[q].parammicro->degradable=1;
                        Inter[q].parammicro->alpha=propinter[j].alpha;
                        Inter[q].parammicro->gamma=propinter[j].gamma;
                        Inter[q].parammicro->Yo=propinter[j].Yo;
                        Inter[q].parammicro->Yc=propinter[j].Yc;
                        Inter[q].parammicro->n=propinter[j].n;
                        Inter[q].parammicro->kn=propinter[j].kn;
                        Inter[q].parammicro->knc=propinter[j].knc;
                        Inter[q].parammicro->kt=propinter[j].kt;

                        Inter[q].comp="Cohesif";
                        break;
                        }
                     }
                  }
                  else{
                     Inter[q].parammicro->degradable=0;
                     }
                  }
               else { //fibre-matrice
                  for(unsigned j=0;j<propinter.size();++j){
                     if (propinter[j].name=="fibre-matrice"){
                        Inter[q].parammicro->Gcrit=propinter[j].Gcrit;
                        Inter[q].coefcontact=propinter[j].coeffrottement;
                        Inter[q].parammicro->degradable=1;
                        Inter[q].parammicro->alpha=propinter[j].alpha;
                        Inter[q].parammicro->gamma=propinter[j].gamma;
                        Inter[q].parammicro->Yo=propinter[j].Yo;
                        Inter[q].parammicro->Yc=propinter[j].Yc;
                        Inter[q].parammicro->n=propinter[j].n;
                        Inter[q].parammicro->kn=propinter[j].kn;
                        Inter[q].parammicro->knc=propinter[j].knc;
                        Inter[q].parammicro->kt=propinter[j].kt;
                        
                        Inter[q].comp="Cohesif";
                        break;
                        }
                     }
                  }
              }
         }
    }
}
