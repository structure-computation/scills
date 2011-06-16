
template<class TV2> void modif_etat_cohesif(TV2 &Inter){
    double eps=1e-3;
    for(unsigned i=0;i<Inter.size();i++){
        if(Inter[i].comp=="Cohesif"){
            double mind=1;
            for(unsigned j=0;j<Inter[i].side[0].nodeeq.size();j++)
                {mind=min(Inter[i].parammicro->d[j],mind);
 //                cout<< "Interface numero "   << i << endl;
//                 cout<< "d[j] "   << Inter[i].parammicro->d[j] <<  " dmax " << //Inter[i].parammicro->dmax[j] <<  endl;
                }


            if(mind<1.-eps){
                Inter[i].parammicro->dmax=Inter[i].parammicro->d;
                cout << "Interface " << i << " passe pas contact" << endl;
                cout<< "Valeur d apres : " <<  Inter[i].parammicro->d <<  endl;

                Inter[i].parammicro->Ymax=Inter[i].parammicro->Y;
                }
            else{
                Inter[i].comp="Contact";
                cout << "Interface " << i << " passe contact" << endl;
                cout<< "Valeur dmax avant : " <<  Inter[i].parammicro->dmax <<  endl;
                Inter[i].parammicro->dmax=Inter[i].parammicro->d;
                cout<< "Valeur dmax apres : " <<  Inter[i].parammicro->dmax <<  endl;
            }
        }
    }
}

