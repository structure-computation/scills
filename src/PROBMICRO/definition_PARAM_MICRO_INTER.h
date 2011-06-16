
#ifndef PARAM_MICRO_INTER_H
#define PARAM_MICRO_INTER_H
using namespace LMT;
using namespace std;

struct PARAM_MICRO_INTER{
    PARAM_MICRO_INTER(unsigned nbnodeeq){Ymax.resize(nbnodeeq,0.);Y.resize(nbnodeeq,0.);d.resize(nbnodeeq,0.);dmax.resize(nbnodeeq,0.);}
    double Gcrit;
    bool degradable;
    Vec<double> dmax,d;
    Vec<double> Y,Ymax;
    double kn, kt,knc;
    double gamma,alpha, Yc, Yo, n;
};

#endif //PARAM_MICRO_INTER_H

