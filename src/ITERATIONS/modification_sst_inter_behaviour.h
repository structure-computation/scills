
void modification_sst_inter_behaviour(Vec<Sst> &S, Vec<Interface> &Inter, TimeParameters &param_incr){
    double eps=1e-3;
    for(unsigned i=0;i<Inter.size();i++){
        if(Inter[i].comp=="Cohesive"){
            double mind=1;
            for(unsigned j=0;j<Inter[i].side[0].nodeeq.size();j++)
                mind=min(Inter[i].param_comp->t[1].d[j],mind);
                std::cout << Inter[i].param_comp->t[1].d << "  "  << Inter[i].param_comp->dmax <<" min " << mind<< endl;

            if(mind<1.-eps){
                Inter[i].param_comp->dmax=Inter[i].param_comp->t[1].d;
                std::cout << "Interface " << i << " passe pas contact" << endl;
                std::cout<< "Valeur d apres : " <<  Inter[i].param_comp->t[1].d <<  endl;

                Inter[i].param_comp->t[1].Ymax=Inter[i].param_comp->t[1].Y;
            }else{
                Inter[i].comp="Contact";
                std::cout << "Interface " << i << " passe contact" << endl;
                std::cout<< "Valeur dmax avant : " <<  Inter[i].param_comp->t[1].dmax <<  endl;
                Inter[i].param_comp->t[1].dmax=Inter[i].param_comp->t[1].d;
                std::cout<< "Valeur dmax apres : " <<  Inter[i].param_comp->t[1].dmax <<  endl;
            }
        }
    }
}

