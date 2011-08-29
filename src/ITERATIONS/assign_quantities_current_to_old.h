#ifndef ASS_QUANTITIES_H
#define ASS_QUANTITIES_H

template<class TV1,class TV2>
void assign_quantities_current_to_old(TV1 &S, TV2 &Inter, Param &process){

//     if (process.latin->save_depl_SST==1)
    PRINT(process.temps->pt_cur);
    
    for(unsigned i=0;i<S.size();i++){
        //PRINT(S[i].t_post.size());
        //PRINT(S[i].t.size());
        //    S[i].t[0].q=S[i].t[1].q;
        //S[i].t[0].oldq=S[i].t[1].oldq;
        S[i].t_post[process.temps->pt_cur]=S[i].t[1];
    }
    PRINT("passe--");
    for(unsigned i=0;i<Inter.size();i++){
        for(unsigned j=0;j<Inter[i].side.size();j++){
            Inter[i].side[j].t[0].F=Inter[i].side[j].t[1].F;
            Inter[i].side[j].t[0].Wp=Inter[i].side[j].t[1].Wp;
            Inter[i].side[j].t[0].Fchap=Inter[i].side[j].t[1].Fchap;
            Inter[i].side[j].t[0].Wpchap=Inter[i].side[j].t[1].Wpchap;
            Inter[i].side[j].t[0].W=Inter[i].side[j].t[1].W;
            Inter[i].side[j].t[0].Wchap=Inter[i].side[j].t[1].Wchap;
            Inter[i].side[j].t_post[process.temps->pt_cur]=Inter[i].side[j].t[1];
      //       if (process.reprise_calcul>0 and process.temps->pt_cur!= (int)process.temps->nbpastemps){//vaut mieux laisser reprendre a partir de la solution en cours car plus proche du resultat...
      //         Inter[i].side[j].t[1]       = Inter[i].side[j].t_post[process.temps->pt_cur+1];
      //         Inter[i].side[j].t[1].oldF  = Inter[i].side[j].t_post[process.temps->pt_cur+1].F;
      //         Inter[i].side[j].t[1].oldWp = Inter[i].side[j].t_post[process.temps->pt_cur+1].Wp;
      //         Inter[i].side[j].t[1].oldW  = Inter[i].side[j].t_post[process.temps->pt_cur+1].W;
      //       }

        }
    }
}

template<class TV2>
void assign_t_post(TV2 &Inter, Param &process){

for(unsigned i=0;i<Inter.size();i++){
   for(unsigned j=0;j<Inter[i].side.size();j++){
      Inter[i].side[j].t_post[process.temps->pt_cur]=Inter[i].side[j].t[1];
   }
}

}
#endif //ASS_QUANTITIES_H
