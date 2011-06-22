
#include "mesh/bar.h"
#include "mesh/triangle.h"
#include "mesh/quad.h"
#include "mesh/bar_3.h"
#include "mesh/triangle_6.h"
#include "mesh/quad_8.h"

namespace LMT {
struct add_elem_Ne{

 /// pour les Bar
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Bar,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=e.pos(1)[0]-e.pos(0)[0]; double reg1=e.pos(1)[1]-e.pos(0)[1]; double reg2=pow(reg1,2); double reg3=pow(reg0,2); reg2=reg3+reg2;
   reg2=pow(reg2,0.5); reg3=reg0/reg2; reg2=reg1/reg2; reg3=reg0*reg3; reg2=reg1*reg2;
   reg2=reg3+reg2; Ne(0,0)+=0.5*reg2; Ne(0,2)+=0.5*reg2; Ne(1,1)+=0.5*reg2; Ne(1,3)+=0.5*reg2;

}
 /// pour les Triangle
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Triangle,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=e.pos(1)[0]-e.pos(0)[0]; double reg1=e.pos(1)[1]-e.pos(0)[1]; double reg2=pow(reg1,2); double reg3=pow(reg0,2); double reg4=e.pos(1)[2]-e.pos(0)[2];
   double reg5=pow(reg4,2); reg2=reg3+reg2; reg5=reg2+reg5; reg5=pow(reg5,0.5); reg2=e.pos(2)[1]-e.pos(0)[1];
   reg3=e.pos(2)[0]-e.pos(0)[0]; double reg6=reg1/reg5; double reg7=reg0/reg5; reg5=reg4/reg5; double reg8=e.pos(2)[2]-e.pos(0)[2];
   double reg9=reg7*reg3; double reg10=reg6*reg2; reg10=reg9+reg10; reg9=reg5*reg8; reg9=reg10+reg9;
   reg10=reg7*reg9; double reg11=reg6*reg9; reg11=reg2-reg11; reg10=reg3-reg10; double reg12=reg5*reg9;
   reg12=reg8-reg12; double reg13=pow(reg10,2); double reg14=pow(reg11,2); reg14=reg13+reg14; reg13=pow(reg12,2);
   reg13=reg14+reg13; reg13=pow(reg13,0.5); reg10=reg10/reg13; reg11=reg11/reg13; reg2=reg2*reg11;
   reg3=reg3*reg10; reg11=reg1*reg11; reg10=reg0*reg10; reg13=reg12/reg13; reg6=reg1*reg6;
   reg7=reg0*reg7; reg8=reg8*reg13; reg2=reg3+reg2; reg13=reg4*reg13; reg11=reg10+reg11;
   reg6=reg7+reg6; reg5=reg4*reg5; reg13=reg11+reg13; reg5=reg6+reg5; reg8=reg2+reg8;
   reg0=reg5*reg8; reg1=reg9*reg13; reg1=reg0-reg1; Ne(0,0)+=0.16666666666666668517*reg1; Ne(1,1)+=0.16666666666666668517*reg1;
   Ne(2,2)+=0.16666666666666668517*reg1; Ne(0,3)+=0.16666666666666665741*reg1; Ne(0,6)+=0.16666666666666665741*reg1; Ne(1,4)+=0.16666666666666665741*reg1; Ne(1,7)+=0.16666666666666665741*reg1;
   Ne(2,5)+=0.16666666666666665741*reg1; Ne(2,8)+=0.16666666666666665741*reg1;

}
 /// pour les Quad
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.21132486540518713447*e.pos(0)[0]; double reg1=0.21132486540518713447*e.pos(1)[0]; double reg2=0.21132486540518713447*e.pos(1)[1]; double reg3=0.21132486540518713447*e.pos(0)[1]; double reg4=0.78867513459481286553*e.pos(2)[1];
   double reg5=0.21132486540518713447*e.pos(0)[2]; double reg6=reg2-reg3; double reg7=0.21132486540518713447*e.pos(1)[2]; double reg8=0.78867513459481286553*e.pos(0)[0]; double reg9=0.78867513459481286553*e.pos(0)[1];
   double reg10=reg1-reg0; double reg11=0.78867513459481286553*e.pos(2)[0]; double reg12=0.78867513459481286553*e.pos(1)[1]; double reg13=0.78867513459481286553*e.pos(1)[0]; double reg14=0.21132486540518713447*e.pos(2)[1];
   double reg15=reg7-reg5; double reg16=0.21132486540518713447*e.pos(2)[0]; double reg17=0.78867513459481286553*e.pos(1)[2]; double reg18=0.78867513459481286553*e.pos(2)[2]; double reg19=0.78867513459481286553*e.pos(0)[2];
   double reg20=0.78867513459481286553*e.pos(3)[1]; reg6=reg4+reg6; double reg21=0.78867513459481286553*e.pos(3)[0]; reg10=reg11+reg10; double reg22=reg12-reg9;
   double reg23=reg13-reg8; double reg24=0.21132486540518713447*e.pos(2)[2]; double reg25=reg17-reg19; reg22=reg22+reg14; reg23=reg16+reg23;
   double reg26=0.21132486540518713447*e.pos(3)[1]; double reg27=0.21132486540518713447*e.pos(3)[0]; reg10=reg10-reg21; reg6=reg6-reg20; reg15=reg18+reg15;
   double reg28=0.78867513459481286553*e.pos(3)[2]; double reg29=pow(reg6,2); reg23=reg23-reg27; reg22=reg22-reg26; reg15=reg15-reg28;
   double reg30=pow(reg10,2); reg25=reg24+reg25; double reg31=0.21132486540518713447*e.pos(3)[2]; double reg32=pow(reg15,2); reg25=reg25-reg31;
   reg29=reg30+reg29; reg30=pow(reg22,2); double reg33=pow(reg23,2); reg9=reg2+reg9; reg2=pow(reg25,2);
   reg32=reg29+reg32; reg8=reg1+reg8; reg12=reg3+reg12; reg30=reg33+reg30; reg13=reg0+reg13;
   reg8=reg16-reg8; reg9=reg14-reg9; reg19=reg7+reg19; reg17=reg5+reg17; reg2=reg30+reg2;
   reg12=reg4-reg12; reg32=pow(reg32,0.5); reg13=reg11-reg13; reg9=reg20+reg9; reg19=reg24-reg19;
   reg8=reg21+reg8; reg2=pow(reg2,0.5); reg0=reg6/reg32; reg1=reg10/reg32; reg27=reg13+reg27;
   reg26=reg12+reg26; reg17=reg18-reg17; reg32=reg15/reg32; reg3=reg1*reg27; reg4=reg0*reg26;
   reg19=reg28+reg19; reg31=reg17+reg31; reg5=reg0*reg9; reg7=reg1*reg8; reg11=reg23/reg2;
   reg12=reg22/reg2; reg4=reg3+reg4; reg3=reg32*reg19; reg13=reg32*reg31; reg5=reg7+reg5;
   reg7=reg27*reg11; reg2=reg25/reg2; reg14=reg26*reg12; reg16=reg9*reg12; reg17=reg8*reg11;
   reg3=reg5+reg3; reg5=reg31*reg2; reg13=reg4+reg13; reg14=reg7+reg14; reg4=reg0*reg3;
   reg7=reg1*reg3; reg5=reg14+reg5; reg16=reg17+reg16; reg14=reg19*reg2; reg17=reg0*reg13;
   reg18=reg1*reg13; reg18=reg27-reg18; reg20=reg12*reg5; reg21=reg32*reg3; reg4=reg9-reg4;
   reg24=reg32*reg13; reg7=reg8-reg7; reg14=reg16+reg14; reg17=reg26-reg17; reg16=reg11*reg5;
   reg28=pow(reg18,2); reg29=reg12*reg14; reg30=reg11*reg14; reg33=reg2*reg5; reg20=reg26-reg20;
   reg16=reg27-reg16; reg21=reg19-reg21; reg24=reg31-reg24; double reg34=pow(reg4,2); double reg35=pow(reg17,2);
   double reg36=pow(reg7,2); double reg37=pow(reg20,2); double reg38=pow(reg16,2); reg33=reg31-reg33; double reg39=reg2*reg14;
   reg30=reg8-reg30; reg29=reg9-reg29; double reg40=pow(reg24,2); reg35=reg28+reg35; reg34=reg36+reg34;
   reg28=pow(reg21,2); reg36=pow(reg29,2); reg40=reg35+reg40; reg39=reg19-reg39; reg35=pow(reg33,2);
   reg28=reg34+reg28; reg37=reg38+reg37; reg34=pow(reg30,2); reg36=reg34+reg36; reg34=pow(reg39,2);
   reg28=pow(reg28,0.5); reg35=reg37+reg35; reg40=pow(reg40,0.5); reg4=reg4/reg28; reg7=reg7/reg28;
   reg17=reg17/reg40; reg35=pow(reg35,0.5); reg34=reg36+reg34; reg18=reg18/reg40; reg40=reg24/reg40;
   reg24=reg26*reg17; reg36=reg10*reg18; reg18=reg27*reg18; reg16=reg16/reg35; reg17=reg6*reg17;
   reg20=reg20/reg35; reg0=reg6*reg0; reg37=reg10*reg7; reg6=reg6*reg4; reg34=pow(reg34,0.5);
   reg28=reg21/reg28; reg4=reg9*reg4; reg7=reg8*reg7; reg1=reg10*reg1; reg12=reg22*reg12;
   reg10=reg22*reg20; reg21=reg23*reg16; reg6=reg37+reg6; reg37=reg15*reg28; reg35=reg33/reg35;
   reg20=reg26*reg20; reg16=reg27*reg16; reg29=reg29/reg34; reg30=reg30/reg34; reg0=reg1+reg0;
   reg28=reg19*reg28; reg1=reg15*reg40; reg24=reg18+reg24; reg4=reg7+reg4; reg40=reg31*reg40;
   reg32=reg15*reg32; reg17=reg36+reg17; reg11=reg23*reg11; reg8=reg8*reg30; reg9=reg9*reg29;
   reg34=reg39/reg34; reg30=reg23*reg30; reg32=reg0+reg32; reg29=reg22*reg29; reg28=reg4+reg28;
   reg0=reg25*reg35; reg10=reg21+reg10; reg37=reg6+reg37; reg12=reg11+reg12; reg35=reg31*reg35;
   reg20=reg16+reg20; reg40=reg24+reg40; reg1=reg17+reg1; reg2=reg25*reg2; reg9=reg8+reg9;
   reg28=reg32*reg28; reg19=reg19*reg34; reg37=reg3*reg37; reg1=reg13*reg1; reg2=reg12+reg2;
   reg40=reg32*reg40; reg35=reg20+reg35; reg29=reg30+reg29; reg34=reg25*reg34; reg0=reg10+reg0;
   reg19=reg9+reg19; reg34=reg29+reg34; reg37=reg28-reg37; reg1=reg40-reg1; reg35=reg2*reg35;
   reg0=reg5*reg0; reg3=0.1555021169820365473*reg37; reg4=0.1555021169820365473*reg1; reg5=0.011164549684630114537*reg37; reg6=0.04166666666666666908*reg1;
   reg34=reg14*reg34; reg0=reg35-reg0; reg19=reg2*reg19; reg2=0.011164549684630114537*reg1; reg7=0.04166666666666666908*reg37;
   reg8=0.04166666666666666908*reg0; reg34=reg19-reg34; reg9=0.011164549684630114537*reg0; reg3=reg6+reg3; reg5=reg6+reg5;
   reg6=0.1555021169820365473*reg0; reg2=reg2+reg7; reg4=reg7+reg4; reg9=reg3+reg9; reg3=0.011164549684630114537*reg34;
   reg4=reg8+reg4; reg7=0.04166666666666666908*reg34; reg6=reg5+reg6; reg5=0.1555021169820365473*reg34; reg8=reg2+reg8;
   Ne(0,3)+=reg6+reg7; Ne(1,4)+=reg6+reg7; Ne(2,5)+=reg6+reg7; Ne(0,6)+=reg4+reg3; Ne(1,7)+=reg4+reg3;
   Ne(2,8)+=reg4+reg3; Ne(0,0)+=reg8+reg5; Ne(1,1)+=reg8+reg5; Ne(2,2)+=reg8+reg5; Ne(0,9)+=reg7+reg9;
   Ne(1,10)+=reg7+reg9; Ne(2,11)+=reg7+reg9;

}
 /// pour les Bar_3
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Bar_3,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=2.1547005383792514621*e.pos(0)[1]; double reg1=0.15470053837925146212*e.pos(1)[1]; double reg2=2.1547005383792514621*e.pos(1)[0]; double reg3=0.15470053837925146212*e.pos(0)[0]; double reg4=2.1547005383792514621*e.pos(1)[1];
   double reg5=2.1547005383792514621*e.pos(0)[0]; double reg6=0.15470053837925146212*e.pos(1)[0]; double reg7=0.15470053837925146212*e.pos(0)[1]; double reg8=2.3094010767585029242*e.pos(2)[1]; reg4=reg4+reg7;
   reg0=reg1+reg0; double reg9=2.3094010767585029242*e.pos(2)[0]; reg5=reg6+reg5; reg2=reg2+reg3; reg4=reg4-reg8;
   reg2=reg2-reg9; double reg10=reg8-reg0; double reg11=reg9-reg5; double reg12=pow(reg4,2); double reg13=pow(reg10,2);
   double reg14=pow(reg2,2); double reg15=pow(reg11,2); reg12=reg14+reg12; reg13=reg15+reg13; reg12=pow(reg12,0.5);
   reg13=pow(reg13,0.5); reg14=reg10/reg13; reg13=reg11/reg13; reg15=reg2/reg12; reg12=reg4/reg12;
   reg14=reg10*reg14; reg13=reg11*reg13; reg15=reg2*reg15; reg12=reg4*reg12; reg12=reg15+reg12;
   reg14=reg13+reg14; reg2=0.22767090063073975644*reg14; reg4=0.061004233964073109089*reg12; reg10=0.061004233964073109089*reg14; reg11=0.22767090063073975644*reg12;
   reg13=0.33333333333333335264*reg14; reg15=0.33333333333333335264*reg12; Ne(0,0)+=reg2-reg4; Ne(1,1)+=reg2-reg4; Ne(0,2)+=reg11-reg10;
   Ne(1,3)+=reg11-reg10; Ne(0,4)+=reg13+reg15; Ne(1,5)+=reg13+reg15;

}
 /// pour les Triangle_6
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Triangle_6,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.78379396366385999995*e.pos(0)[1]; double reg1=0.56758792732771999991*e.pos(1)[1]; double reg2=0.78379396366385999995*e.pos(0)[0]; double reg3=0.56758792732771999991*e.pos(1)[0]; double reg4=0.56758792732771999991*e.pos(0)[1];
   double reg5=0.78379396366385999995*e.pos(1)[1]; double reg6=0.56758792732771999991*e.pos(0)[0]; double reg7=0.78379396366385999995*e.pos(1)[0]; double reg8=reg5+reg4; double reg9=1.3513818909915799999*e.pos(3)[1];
   double reg10=reg7+reg6; double reg11=1.3513818909915799999*e.pos(3)[0]; double reg12=0.78379396366385999995*e.pos(0)[2]; double reg13=0.56758792732771999991*e.pos(1)[2]; double reg14=0.78379396366385999995*e.pos(1)[2];
   double reg15=0.56758792732771999991*e.pos(0)[2]; double reg16=reg3+reg2; double reg17=reg1+reg0; double reg18=reg14+reg15; double reg19=1.3513818909915799999*e.pos(3)[2];
   double reg20=reg13+reg12; double reg21=reg11-reg16; double reg22=reg9-reg17; double reg23=1.78379396366386*e.pos(4)[1]; reg8=reg8-reg9;
   double reg24=0.63369514596091600003*e.pos(1)[0]; double reg25=1.78379396366386*e.pos(4)[0]; reg10=reg10-reg11; double reg26=2.2673902919218320001*e.pos(0)[0]; double reg27=0.63369514596091600003*e.pos(1)[1];
   double reg28=2.2673902919218320001*e.pos(0)[1]; double reg29=0.63369514596091600003*e.pos(0)[0]; double reg30=2.2673902919218320001*e.pos(1)[0]; reg21=reg25+reg21; double reg31=0.63369514596091600003*e.pos(1)[2];
   double reg32=2.9010854378827480001*e.pos(3)[1]; double reg33=reg27+reg28; double reg34=reg24+reg26; double reg35=2.2673902919218320001*e.pos(1)[1]; double reg36=2.9010854378827480001*e.pos(3)[0];
   double reg37=0.63369514596091600003*e.pos(0)[1]; double reg38=2.2673902919218320001*e.pos(0)[2]; double reg39=0.43241207267228000009*e.pos(4)[1]; reg5=reg5-reg0; double reg40=0.43241207267228000009*e.pos(4)[0];
   reg10=reg10+reg25; double reg41=1.78379396366386*e.pos(5)[0]; double reg42=reg19-reg20; reg8=reg8+reg23; double reg43=1.78379396366386*e.pos(5)[1];
   reg18=reg18-reg19; double reg44=1.78379396366386*e.pos(4)[2]; reg7=reg7-reg2; reg22=reg23+reg22; double reg45=2.9010854378827480001*e.pos(3)[2];
   double reg46=reg31+reg38; reg7=reg40+reg7; double reg47=0.43241207267228000009*e.pos(5)[0]; double reg48=0.43241207267228000009*e.pos(4)[2]; double reg49=0.366304854039084*e.pos(4)[1];
   double reg50=reg32-reg33; reg5=reg39+reg5; double reg51=0.43241207267228000009*e.pos(5)[1]; double reg52=0.366304854039084*e.pos(4)[0]; double reg53=reg36-reg34;
   reg14=reg14-reg12; reg42=reg44+reg42; reg22=reg22-reg43; reg21=reg21-reg41; reg41=reg10-reg41;
   reg43=reg8-reg43; reg18=reg18+reg44; reg8=1.78379396366386*e.pos(5)[2]; reg10=0.63369514596091600003*e.pos(0)[2]; double reg54=2.2673902919218320001*e.pos(1)[2];
   reg35=reg35+reg37; reg30=reg30+reg29; double reg55=0.78379396366385999995*e.pos(2)[0]; reg47=reg7-reg47; reg7=2.710505431213761085e-20*e.pos(3)[0];
   double reg56=reg29-reg24; double reg57=pow(reg21,2); reg18=reg18-reg8; double reg58=pow(reg22,2); double reg59=pow(reg43,2);
   double reg60=pow(reg41,2); reg51=reg5-reg51; reg8=reg42-reg8; reg5=reg37-reg27; reg42=2.710505431213761085e-20*e.pos(3)[1];
   double reg61=0.43241207267228000009*e.pos(5)[2]; reg14=reg48+reg14; reg35=reg35-reg32; double reg62=0.366304854039084*e.pos(5)[1]; reg50=reg50+reg49;
   double reg63=0.78379396366385999995*e.pos(2)[1]; reg30=reg30-reg36; double reg64=reg45-reg46; double reg65=0.366304854039084*e.pos(5)[0]; reg54=reg54+reg10;
   reg53=reg53+reg52; double reg66=0.366304854039084*e.pos(4)[2]; double reg67=0.78379396366385999995*e.pos(2)[2]; reg61=reg14-reg61; reg5=reg5-reg42;
   reg14=0.56758792732771999991*e.pos(2)[0]; double reg68=0.56758792732771999991*e.pos(2)[1]; double reg69=pow(reg51,2); reg50=reg50-reg62; reg54=reg54-reg45;
   double reg70=pow(reg8,2); double reg71=1.78379396366386*e.pos(3)[0]; reg35=reg49+reg35; double reg72=1.78379396366386*e.pos(3)[1]; double reg73=reg4+reg63;
   double reg74=3.2673902919218320001*e.pos(4)[1]; double reg75=pow(reg47,2); double reg76=3.2673902919218320001*e.pos(4)[0]; reg64=reg64+reg66; reg53=reg53-reg65;
   reg56=reg56-reg7; reg30=reg52+reg30; double reg77=2.710505431213761085e-20*e.pos(3)[2]; double reg78=pow(reg18,2); double reg79=0.366304854039084*e.pos(5)[2];
   double reg80=reg10-reg31; double reg81=reg6+reg55; reg58=reg57+reg58; reg59=reg60+reg59; reg80=reg80-reg77;
   reg69=reg75+reg69; reg5=reg74+reg5; reg57=3.2673902919218320001*e.pos(5)[1]; reg60=pow(reg61,2); reg75=reg2+reg14;
   double reg82=reg0+reg68; double reg83=0.56758792732771999991*e.pos(2)[2]; double reg84=0.63369514596091600003*e.pos(2)[1]; double reg85=0.63369514596091600003*e.pos(2)[0]; reg64=reg64-reg79;
   double reg86=3.2673902919218320001*e.pos(5)[0]; reg56=reg76+reg56; double reg87=3.2673902919218320001*e.pos(4)[2]; reg54=reg66+reg54; double reg88=pow(reg53,2);
   reg62=reg35-reg62; reg35=pow(reg50,2); reg65=reg30-reg65; reg30=1.78379396366386*e.pos(3)[2]; double reg89=reg15+reg67;
   double reg90=0.43241207267228000009*e.pos(3)[1]; reg0=reg63-reg0; reg78=reg59+reg78; reg59=0.43241207267228000009*e.pos(3)[0]; reg2=reg55-reg2;
   reg81=reg81-reg71; reg73=reg73-reg72; reg70=reg58+reg70; reg89=reg89-reg30; reg55=0.43241207267228000009*e.pos(3)[2];
   reg67=reg67-reg12; reg73=reg23+reg73; reg72=reg72+reg82; reg90=reg0-reg90; reg0=pow(reg62,2);
   reg58=pow(reg65,2); reg63=1.3513818909915799999*e.pos(5)[1]; reg59=reg2-reg59; reg71=reg71+reg75; reg2=1.3513818909915799999*e.pos(5)[0];
   reg79=reg54-reg79; reg60=reg69+reg60; reg81=reg25+reg81; reg54=0.366304854039084*e.pos(3)[0]; reg57=reg5-reg57;
   reg28=reg28+reg84; reg86=reg56-reg86; reg78=pow(reg78,0.5); reg80=reg87+reg80; reg5=0.366304854039084*e.pos(3)[1];
   reg56=3.2673902919218320001*e.pos(5)[2]; reg69=0.63369514596091600003*e.pos(2)[2]; reg26=reg26+reg85; double reg91=pow(reg64,2); reg35=reg88+reg35;
   reg70=pow(reg70,0.5); reg12=reg12+reg83; reg91=reg35+reg91; reg25=reg25-reg71; reg35=reg21/reg70;
   reg88=reg26+reg54; reg73=reg73-reg63; double reg92=reg22/reg70; double reg93=reg28+reg5; reg38=reg38+reg69;
   reg23=reg23-reg72; double reg94=0.366304854039084*e.pos(3)[2]; double reg95=1.3513818909915799999*e.pos(5)[2]; reg30=reg30+reg12; reg89=reg44+reg89;
   double reg96=pow(reg86,2); double reg97=reg41/reg78; double reg98=reg43/reg78; double reg99=pow(reg57,2); reg56=reg80-reg56;
   reg80=3.2673902919218320001*e.pos(3)[1]; double reg100=reg37-reg84; double reg101=2.2673902919218320001*e.pos(2)[0]; reg55=reg67-reg55; reg67=3.2673902919218320001*e.pos(3)[0];
   double reg102=reg29-reg85; reg39=reg90+reg39; reg90=2.2673902919218320001*e.pos(2)[1]; reg40=reg59+reg40; reg59=pow(reg79,2);
   reg60=pow(reg60,0.5); reg81=reg81-reg2; reg0=reg58+reg0; reg101=reg29+reg101; reg90=reg37+reg90;
   reg58=pow(reg56,2); reg91=pow(reg91,0.5); reg99=reg96+reg99; reg96=reg52-reg88; double reg103=2.9010854378827480001*e.pos(5)[0];
   double reg104=2.2673902919218320001*e.pos(2)[2]; reg89=reg89-reg95; double reg105=reg49-reg93; double reg106=2.9010854378827480001*e.pos(5)[1]; double reg107=3.2673902919218320001*e.pos(3)[2];
   double reg108=reg10-reg69; reg78=reg18/reg78; double reg109=reg38+reg94; reg80=reg100-reg80; reg67=reg102-reg67;
   reg59=reg0+reg59; reg0=reg98*reg73; reg100=reg97*reg81; reg44=reg44-reg30; reg23=reg63+reg23;
   reg48=reg55+reg48; reg25=reg2+reg25; reg70=reg8/reg70; reg55=reg92*reg39; reg102=reg51/reg60;
   double reg110=reg35*reg40; double reg111=reg47/reg60; double reg112=5.42101086242752217e-20*e.pos(5)[1]; reg74=reg80+reg74; reg80=reg70*reg48;
   reg5=reg90-reg5; reg90=2.9010854378827480001*e.pos(5)[2]; double reg113=reg66-reg109; reg44=reg95+reg44; double reg114=5.42101086242752217e-20*e.pos(5)[0];
   double reg115=reg53/reg91; double reg116=reg50/reg91; reg104=reg10+reg104; reg54=reg101-reg54; reg101=reg78*reg89;
   reg76=reg67+reg76; reg96=reg96+reg103; reg105=reg105+reg106; reg67=reg102*reg23; double reg117=reg111*reg25;
   reg55=reg110+reg55; reg58=reg99+reg58; reg59=pow(reg59,0.5); reg107=reg108-reg107; reg60=reg61/reg60;
   reg0=reg100+reg0; reg91=reg64/reg91; reg99=reg62/reg59; reg80=reg55+reg80; reg55=reg60*reg44;
   reg54=reg52+reg54; reg76=reg76-reg114; reg101=reg0+reg101; reg94=reg104-reg94; reg0=reg116*reg105;
   reg52=5.42101086242752217e-20*e.pos(5)[2]; reg5=reg49+reg5; reg113=reg113+reg90; reg49=reg115*reg96; reg87=reg107+reg87;
   reg100=reg65/reg59; reg58=pow(reg58,0.5); reg67=reg117+reg67; reg74=reg74-reg112; reg94=reg66+reg94;
   reg59=reg79/reg59; reg0=reg49+reg0; reg49=reg35*reg80; reg55=reg67+reg55; reg87=reg87-reg52;
   reg66=reg97*reg101; reg67=reg92*reg80; reg104=reg98*reg101; reg5=reg5-reg106; reg107=reg100*reg76;
   reg108=reg91*reg113; reg54=reg54-reg103; reg110=reg57/reg58; reg117=reg86/reg58; double reg118=reg99*reg74;
   double reg119=reg111*reg55; reg108=reg0+reg108; reg49=reg40-reg49; reg0=reg102*reg55; reg94=reg94-reg90;
   reg67=reg39-reg67; reg58=reg56/reg58; reg66=reg81-reg66; double reg120=reg110*reg5; reg104=reg73-reg104;
   double reg121=reg78*reg101; reg118=reg107+reg118; reg107=reg59*reg87; double reg122=reg117*reg54; double reg123=reg70*reg80;
   double reg124=reg58*reg94; reg123=reg48-reg123; reg120=reg122+reg120; reg122=reg116*reg108; double reg125=pow(reg49,2);
   reg121=reg89-reg121; double reg126=pow(reg104,2); double reg127=pow(reg66,2); double reg128=pow(reg67,2); double reg129=reg115*reg108;
   reg119=reg25-reg119; reg107=reg118+reg107; reg0=reg23-reg0; reg118=reg60*reg55; double reg130=pow(reg121,2);
   reg126=reg127+reg126; reg125=reg128+reg125; reg124=reg120+reg124; reg118=reg44-reg118; reg120=reg100*reg107;
   reg127=pow(reg0,2); reg128=pow(reg119,2); double reg131=pow(reg123,2); reg129=reg96-reg129; double reg132=reg99*reg107;
   double reg133=reg91*reg108; reg122=reg105-reg122; double reg134=reg110*reg124; reg125=reg131+reg125; reg120=reg76-reg120;
   reg132=reg74-reg132; reg131=reg59*reg107; reg130=reg126+reg130; reg126=pow(reg129,2); double reg135=reg117*reg124;
   double reg136=pow(reg122,2); double reg137=pow(reg118,2); reg127=reg128+reg127; reg133=reg113-reg133; reg125=pow(reg125,0.5);
   reg131=reg87-reg131; reg128=pow(reg132,2); reg136=reg126+reg136; reg126=pow(reg120,2); double reg138=reg58*reg124;
   double reg139=pow(reg133,2); reg134=reg5-reg134; reg137=reg127+reg137; reg135=reg54-reg135; reg130=pow(reg130,0.5);
   reg128=reg126+reg128; reg126=pow(reg131,2); reg127=pow(reg134,2); reg49=reg49/reg125; reg137=pow(reg137,0.5);
   reg66=reg66/reg130; double reg140=pow(reg135,2); reg139=reg136+reg139; reg67=reg67/reg125; reg138=reg94-reg138;
   reg104=reg104/reg130; reg0=reg0/reg137; reg136=pow(reg138,2); reg81=reg81*reg66; reg73=reg73*reg104;
   reg130=reg121/reg130; reg125=reg123/reg125; reg66=reg41*reg66; reg98=reg43*reg98; reg104=reg43*reg104;
   reg127=reg140+reg127; reg119=reg119/reg137; reg43=reg22*reg67; reg139=pow(reg139,0.5); reg97=reg41*reg97;
   reg41=reg21*reg49; reg67=reg39*reg67; reg35=reg21*reg35; reg92=reg22*reg92; reg49=reg40*reg49;
   reg126=reg128+reg126; reg136=reg127+reg136; reg21=reg51*reg0; reg22=reg47*reg119; reg137=reg118/reg137;
   reg126=pow(reg126,0.5); reg48=reg48*reg125; reg70=reg8*reg70; reg111=reg47*reg111; reg102=reg51*reg102;
   reg92=reg35+reg92; reg67=reg49+reg67; reg122=reg122/reg139; reg129=reg129/reg139; reg78=reg18*reg78;
   reg119=reg25*reg119; reg0=reg23*reg0; reg73=reg81+reg73; reg89=reg89*reg130; reg125=reg8*reg125;
   reg43=reg41+reg43; reg98=reg97+reg98; reg104=reg66+reg104; reg130=reg18*reg130; reg48=reg67+reg48;
   reg21=reg22+reg21; reg8=reg61*reg137; reg136=pow(reg136,0.5); reg137=reg44*reg137; reg0=reg119+reg0;
   reg130=reg104+reg130; reg89=reg73+reg89; reg115=reg53*reg115; reg116=reg50*reg116; reg96=reg96*reg129;
   reg105=reg105*reg122; reg139=reg133/reg139; reg129=reg53*reg129; reg122=reg50*reg122; reg60=reg61*reg60;
   reg102=reg111+reg102; reg78=reg98+reg78; reg70=reg92+reg70; reg125=reg43+reg125; reg120=reg120/reg126;
   reg132=reg132/reg126; reg135=reg135/reg136; reg134=reg134/reg136; reg113=reg113*reg139; reg8=reg21+reg8;
   reg122=reg129+reg122; reg139=reg64*reg139; reg137=reg0+reg137; reg105=reg96+reg105; reg60=reg102+reg60;
   reg100=reg65*reg100; reg99=reg62*reg99; reg91=reg64*reg91; reg116=reg115+reg116; reg125=reg80*reg125;
   reg76=reg76*reg120; reg74=reg74*reg132; reg48=reg70*reg48; reg132=reg62*reg132; reg120=reg65*reg120;
   reg130=reg101*reg130; reg126=reg131/reg126; reg89=reg78*reg89; reg113=reg105+reg113; reg8=reg55*reg8;
   reg139=reg122+reg139; reg137=reg60*reg137; reg54=reg54*reg135; reg5=reg5*reg134; reg130=reg89-reg130;
   reg0=reg79*reg126; reg132=reg120+reg132; reg99=reg100+reg99; reg59=reg79*reg59; reg136=reg138/reg136;
   reg91=reg116+reg91; reg135=reg86*reg135; reg134=reg57*reg134; reg74=reg76+reg74; reg126=reg87*reg126;
   reg110=reg57*reg110; reg117=reg86*reg117; reg125=reg48-reg125; reg8=reg137-reg8; reg5=reg54+reg5;
   reg94=reg94*reg136; reg134=reg135+reg134; reg136=reg56*reg136; reg18=0.005384432036113586778*reg130; reg21=0.009463616120767210603*reg125;
   reg22=0.021537728144454347112*reg130; reg23=0.021537728144454347112*reg125; reg25=0.088847818743090689935*reg130; reg35=0.088847818743090689935*reg125; reg126=reg74+reg126;
   reg39=0.009463616120767210603*reg130; reg113=reg91*reg113; reg139=reg108*reg139; reg0=reg132+reg0; reg110=reg117+reg110;
   reg59=reg99+reg59; reg40=0.005384432036113586778*reg125; reg58=reg56*reg58; reg126=reg59*reg126; reg21=reg18+reg21;
   reg35=reg22+reg35; reg18=reg40+reg18; reg41=0.009463616120767210603*reg8; reg43=0.021537728144454347112*reg8; reg25=reg23+reg25;
   reg0=reg107*reg0; reg23=reg22+reg23; reg22=0.088847818743090689935*reg8; reg44=0.005384432036113586778*reg8; reg136=reg134+reg136;
   reg139=reg113-reg139; reg58=reg110+reg58; reg94=reg5+reg94; reg40=reg39+reg40; reg5=0.016449618187943419918*reg139;
   reg22=reg23+reg22; reg40=reg40+reg44; reg0=reg126-reg0; reg23=0.028457289286966203713*reg139; reg25=reg25+reg43;
   reg41=reg18+reg41; reg94=reg58*reg94; reg18=0.0018441552587796664112*reg139; reg39=0.0041124045469858549794*reg139; reg136=reg124*reg136;
   reg21=reg44+reg21; reg35=reg43+reg35; reg43=0.0018441552587796664109*reg0; reg44=0.004112404546985854979*reg0; reg18=reg25+reg18;
   reg25=0.016449618187943419918*reg0; reg35=reg5+reg35; reg40=reg23-reg40; reg136=reg94-reg136; reg23=0.028457289286966203713*reg0;
   reg21=reg21+reg39; reg41=reg39+reg41; reg39=0.0041124045469858549794*reg0; reg47=0.016449618187943419916*reg0; reg5=reg22+reg5;
   reg39=reg41+reg39; reg22=0.028457289286966203713*reg136; reg43=reg35+reg43; reg21=reg23-reg21; reg23=0.0018441552587796664111*reg136;
   reg25=reg18+reg25; reg18=0.016449618187943419918*reg136; reg47=reg5+reg47; reg44=reg40-reg44; reg5=0.0041124045469858549794*reg136;
   Ne(0,15)+=reg18+reg43; Ne(1,16)+=reg18+reg43; Ne(2,17)+=reg18+reg43; Ne(0,0)+=reg44-reg5; Ne(1,1)+=reg44-reg5;
   Ne(2,2)+=reg44-reg5; Ne(0,12)+=reg25+reg18; Ne(1,13)+=reg25+reg18; Ne(2,14)+=reg25+reg18; Ne(0,3)+=reg21-reg5;
   Ne(1,4)+=reg21-reg5; Ne(2,5)+=reg21-reg5; Ne(0,9)+=reg47+reg23; Ne(1,10)+=reg47+reg23; Ne(2,11)+=reg47+reg23;
   Ne(0,6)+=reg22-reg39; Ne(1,7)+=reg22-reg39; Ne(2,8)+=reg22-reg39;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(0)[1]; double reg1=0.12200846792814621817*e.pos(1)[1]; double reg2=0.36602540378443865451*e.pos(1)[0]; double reg3=0.12200846792814621817*e.pos(0)[0]; double reg4=0.36602540378443865451*e.pos(1)[1];
   double reg5=0.12200846792814621817*e.pos(0)[1]; double reg6=0.36602540378443865451*e.pos(0)[0]; double reg7=0.12200846792814621817*e.pos(1)[0]; double reg8=1.3660254037844385386*e.pos(2)[0]; double reg9=reg7+reg6;
   double reg10=reg1+reg0; double reg11=1.3660254037844385386*e.pos(2)[1]; double reg12=0.12200846792814621817*e.pos(1)[2]; double reg13=0.36602540378443865451*e.pos(0)[2]; double reg14=reg4+reg5;
   double reg15=0.45534180126147951289*e.pos(2)[1]; double reg16=1.3660254037844385386*e.pos(1)[1]; double reg17=1.3660254037844385386*e.pos(1)[0]; double reg18=0.36602540378443865451*e.pos(1)[2]; double reg19=0.12200846792814621817*e.pos(0)[2];
   double reg20=0.45534180126147951289*e.pos(0)[1]; double reg21=0.45534180126147951289*e.pos(0)[0]; double reg22=0.45534180126147951289*e.pos(2)[0]; double reg23=reg2+reg3; double reg24=reg12+reg13;
   double reg25=1.3660254037844385386*e.pos(2)[2]; double reg26=1.3660254037844385386*e.pos(1)[2]; double reg27=reg14+reg15; double reg28=reg20+reg16; double reg29=0.45534180126147951289*e.pos(1)[0];
   double reg30=1.3660254037844385386*e.pos(3)[1]; double reg31=reg21+reg17; double reg32=0.45534180126147951289*e.pos(1)[1]; double reg33=1.3660254037844385386*e.pos(0)[1]; double reg34=reg18+reg19;
   double reg35=0.45534180126147951289*e.pos(2)[2]; double reg36=0.45534180126147951289*e.pos(0)[2]; double reg37=0.12200846792814621817*e.pos(2)[1]; double reg38=0.12200846792814621817*e.pos(2)[0]; double reg39=reg23+reg22;
   double reg40=0.45534180126147951289*e.pos(3)[0]; double reg41=1.3660254037844385386*e.pos(3)[0]; reg10=reg10+reg11; reg9=reg9+reg8; double reg42=1.3660254037844385386*e.pos(0)[0];
   double reg43=0.45534180126147951289*e.pos(3)[1]; double reg44=0.48803387171258487271*e.pos(4)[0]; reg28=reg37+reg28; double reg45=0.36602540378443865451*e.pos(3)[1]; double reg46=reg27+reg30;
   double reg47=0.36602540378443865451*e.pos(2)[0]; reg31=reg38+reg31; double reg48=reg32+reg33; double reg49=0.36602540378443865451*e.pos(2)[1]; double reg50=0.36602540378443865451*e.pos(3)[0];
   double reg51=reg36+reg26; reg9=reg9+reg40; reg10=reg10+reg43; double reg52=0.48803387171258487271*e.pos(4)[1]; double reg53=reg29+reg42;
   double reg54=0.12200846792814621817*e.pos(2)[2]; double reg55=1.3660254037844385386*e.pos(3)[2]; double reg56=reg34+reg35; reg24=reg24+reg25; double reg57=0.45534180126147951289*e.pos(3)[2];
   double reg58=0.45534180126147951289*e.pos(1)[2]; double reg59=1.3660254037844385386*e.pos(0)[2]; double reg60=reg39+reg41; double reg61=0.12200846792814621817*e.pos(3)[1]; double reg62=0.36602540378443865451*e.pos(2)[2];
   double reg63=reg56+reg55; double reg64=reg53+reg47; double reg65=reg58+reg59; double reg66=reg48+reg49; double reg67=0.66666666666666670528*e.pos(5)[0];
   reg9=reg9-reg44; reg31=reg31+reg50; reg10=reg10-reg52; double reg68=0.66666666666666670528*e.pos(5)[1]; double reg69=reg44-reg60;
   reg24=reg24+reg57; double reg70=0.48803387171258487271*e.pos(4)[2]; double reg71=0.36602540378443865451*e.pos(3)[2]; reg51=reg54+reg51; double reg72=1.8213672050459180516*e.pos(4)[1];
   double reg73=1.8213672050459180516*e.pos(4)[0]; double reg74=reg52-reg46; double reg75=0.12200846792814621817*e.pos(3)[0]; reg28=reg28+reg45; reg31=reg31-reg73;
   reg10=reg10+reg68; double reg76=1.8213672050459180516*e.pos(6)[1]; double reg77=reg70-reg63; reg74=reg68+reg74; reg69=reg67+reg69;
   double reg78=1.8213672050459180516*e.pos(4)[2]; reg51=reg51+reg71; double reg79=reg65+reg62; double reg80=0.12200846792814621817*e.pos(3)[2]; double reg81=reg75+reg64;
   double reg82=1.8213672050459180516*e.pos(6)[0]; reg24=reg24-reg70; double reg83=0.66666666666666670528*e.pos(5)[2]; reg9=reg67+reg9; double reg84=reg61+reg66;
   reg28=reg28-reg72; reg69=reg82+reg69; reg77=reg83+reg77; reg74=reg76+reg74; reg9=reg9-reg82;
   double reg85=0.66666666666666670528*e.pos(7)[0]; reg10=reg10-reg76; double reg86=0.66666666666666670528*e.pos(7)[1]; reg51=reg51-reg78; reg24=reg24+reg83;
   double reg87=1.8213672050459180516*e.pos(6)[2]; double reg88=0.48803387171258487271*e.pos(6)[1]; reg28=reg68+reg28; double reg89=reg6+reg29; double reg90=0.48803387171258487271*e.pos(6)[0];
   reg31=reg67+reg31; double reg91=reg0+reg32; double reg92=reg73-reg81; double reg93=reg2+reg21; double reg94=reg72-reg84;
   double reg95=reg80+reg79; double reg96=reg4+reg20; reg16=reg5+reg16; reg17=reg3+reg17; reg38=reg38+reg93;
   reg77=reg87+reg77; reg3=0.48803387171258487271*e.pos(6)[2]; reg51=reg83+reg51; reg37=reg37+reg96; reg5=reg18+reg36;
   reg69=reg69-reg85; double reg97=reg78-reg95; reg28=reg28-reg88; reg94=reg68+reg94; reg92=reg67+reg92;
   reg31=reg31-reg90; reg74=reg74-reg86; reg11=reg11+reg91; reg9=reg9-reg85; reg67=reg13+reg58;
   reg10=reg10-reg86; reg68=0.66666666666666670528*e.pos(7)[2]; reg8=reg8+reg89; reg24=reg24-reg87; reg28=reg28-reg86;
   reg54=reg54+reg5; reg17=reg22+reg17; reg31=reg31-reg85; reg77=reg77-reg68; reg75=reg8+reg75;
   reg8=0.66666666666666670528*e.pos(4)[0]; reg38=reg41+reg38; reg97=reg83+reg97; reg16=reg15+reg16; reg61=reg11+reg61;
   reg11=0.66666666666666670528*e.pos(4)[1]; reg92=reg90+reg92; reg26=reg19+reg26; reg94=reg88+reg94; reg19=pow(reg74,2);
   reg25=reg25+reg67; reg51=reg51-reg3; reg24=reg24-reg68; reg37=reg30+reg37; reg33=reg1+reg33;
   reg42=reg7+reg42; reg1=pow(reg10,2); reg7=pow(reg69,2); reg30=pow(reg9,2); reg80=reg25+reg80;
   reg25=reg45+reg16; reg41=0.48803387171258487271*e.pos(5)[1]; reg38=reg38-reg8; reg37=reg37-reg11; reg83=0.48803387171258487271*e.pos(5)[0];
   reg85=reg92-reg85; reg86=reg94-reg86; reg92=1.8213672050459180516*e.pos(5)[1]; reg94=0.66666666666666670528*e.pos(4)[2]; reg51=reg51-reg68;
   reg42=reg47+reg42; reg1=reg30+reg1; reg33=reg49+reg33; reg30=pow(reg24,2); reg54=reg55+reg54;
   reg19=reg7+reg19; reg26=reg35+reg26; reg7=pow(reg28,2); reg59=reg12+reg59; reg12=pow(reg31,2);
   reg55=pow(reg77,2); reg75=reg75-reg8; double reg98=1.8213672050459180516*e.pos(5)[0]; double reg99=reg50+reg17; reg97=reg3+reg97;
   reg61=reg61-reg11; reg80=reg80-reg94; double reg100=reg8+reg99; double reg101=0.48803387171258487271*e.pos(5)[2]; double reg102=0.66666666666666670528*e.pos(6)[1];
   reg61=reg61-reg92; double reg103=pow(reg85,2); double reg104=0.66666666666666670528*e.pos(6)[0]; reg75=reg75-reg98; double reg105=pow(reg86,2);
   reg30=reg1+reg30; reg7=reg12+reg7; reg1=pow(reg51,2); reg68=reg97-reg68; reg12=reg11+reg25;
   reg38=reg38-reg83; reg55=reg19+reg55; reg59=reg62+reg59; reg37=reg37-reg41; reg19=reg43+reg33;
   reg97=reg71+reg26; reg54=reg54-reg94; double reg106=reg40+reg42; double reg107=1.8213672050459180516*e.pos(5)[2]; double reg108=pow(reg68,2);
   reg54=reg54-reg101; reg38=reg104+reg38; double reg109=1.8213672050459180516*e.pos(7)[0]; double reg110=reg57+reg59; double reg111=reg92-reg12;
   reg61=reg61+reg102; double reg112=0.48803387171258487271*e.pos(7)[1]; reg80=reg80-reg107; reg1=reg7+reg1; reg7=0.48803387171258487271*e.pos(7)[0];
   reg30=pow(reg30,0.5); reg103=reg105+reg103; reg8=reg8+reg106; reg75=reg75+reg104; reg55=pow(reg55,0.5);
   reg37=reg102+reg37; reg105=1.8213672050459180516*e.pos(7)[1]; reg11=reg11+reg19; double reg113=reg94+reg97; double reg114=reg98-reg100;
   double reg115=0.66666666666666670528*e.pos(6)[2]; double reg116=reg69/reg55; reg111=reg102+reg111; reg94=reg94+reg110; reg1=pow(reg1,0.5);
   double reg117=reg10/reg30; double reg118=reg9/reg30; double reg119=reg41-reg11; double reg120=reg83-reg8; double reg121=reg107-reg113;
   double reg122=reg74/reg55; reg108=reg103+reg108; reg75=reg75-reg7; reg54=reg115+reg54; reg103=1.8213672050459180516*e.pos(7)[2];
   double reg123=0.48803387171258487271*e.pos(7)[2]; reg80=reg80+reg115; reg38=reg38-reg109; reg37=reg37-reg105; reg114=reg104+reg114;
   reg61=reg61-reg112; reg119=reg102+reg119; reg102=reg122*reg37; reg120=reg104+reg120; reg114=reg7+reg114;
   reg54=reg54-reg103; reg121=reg115+reg121; reg111=reg112+reg111; reg80=reg80-reg123; reg104=reg117*reg61;
   double reg124=reg116*reg38; double reg125=reg31/reg1; reg30=reg24/reg30; double reg126=reg28/reg1; double reg127=reg101-reg94;
   reg55=reg77/reg55; double reg128=reg118*reg75; reg108=pow(reg108,0.5); reg1=reg51/reg1; double reg129=reg126*reg111;
   double reg130=reg125*reg114; reg121=reg123+reg121; double reg131=reg30*reg80; reg102=reg124+reg102; reg124=reg55*reg54;
   reg127=reg115+reg127; reg115=reg85/reg108; reg104=reg128+reg104; reg128=reg86/reg108; reg119=reg105+reg119;
   reg120=reg109+reg120; reg131=reg104+reg131; reg124=reg102+reg124; reg102=reg128*reg119; reg104=reg115*reg120;
   reg108=reg68/reg108; double reg132=reg1*reg121; reg127=reg103+reg127; reg129=reg130+reg129; reg130=reg118*reg131;
   double reg133=reg117*reg131; reg102=reg104+reg102; reg104=reg116*reg124; double reg134=reg122*reg124; double reg135=reg108*reg127;
   reg132=reg129+reg132; reg129=reg125*reg132; double reg136=reg30*reg131; reg133=reg61-reg133; reg130=reg75-reg130;
   reg104=reg38-reg104; reg134=reg37-reg134; double reg137=reg55*reg124; double reg138=reg126*reg132; reg135=reg102+reg135;
   reg138=reg111-reg138; reg102=reg1*reg132; double reg139=reg128*reg135; double reg140=reg115*reg135; reg137=reg54-reg137;
   double reg141=pow(reg134,2); double reg142=pow(reg104,2); reg136=reg80-reg136; double reg143=pow(reg133,2); double reg144=pow(reg130,2);
   reg129=reg114-reg129; reg102=reg121-reg102; double reg145=pow(reg136,2); reg143=reg144+reg143; reg144=pow(reg137,2);
   double reg146=pow(reg138,2); double reg147=reg108*reg135; reg139=reg119-reg139; double reg148=pow(reg129,2); reg140=reg120-reg140;
   reg141=reg142+reg141; reg145=reg143+reg145; reg146=reg148+reg146; reg142=pow(reg140,2); reg143=pow(reg139,2);
   reg147=reg127-reg147; reg144=reg141+reg144; reg141=pow(reg102,2); reg145=pow(reg145,0.5); reg144=pow(reg144,0.5);
   reg148=pow(reg147,2); reg141=reg146+reg141; reg143=reg142+reg143; reg134=reg134/reg144; reg104=reg104/reg144;
   reg141=pow(reg141,0.5); reg133=reg133/reg145; reg130=reg130/reg145; reg148=reg143+reg148; reg142=reg9*reg130;
   reg143=reg10*reg133; reg116=reg69*reg116; reg122=reg74*reg122; reg148=pow(reg148,0.5); reg117=reg10*reg117;
   reg118=reg9*reg118; reg145=reg136/reg145; reg133=reg61*reg133; reg130=reg75*reg130; reg38=reg38*reg104;
   reg138=reg138/reg141; reg129=reg129/reg141; reg74=reg74*reg134; reg104=reg69*reg104; reg144=reg137/reg144;
   reg134=reg37*reg134; reg30=reg24*reg30; reg117=reg118+reg117; reg143=reg142+reg143; reg24=reg24*reg145;
   reg141=reg102/reg141; reg145=reg80*reg145; reg111=reg111*reg138; reg114=reg114*reg129; reg134=reg38+reg134;
   reg129=reg31*reg129; reg133=reg130+reg133; reg140=reg140/reg148; reg125=reg31*reg125; reg139=reg139/reg148;
   reg138=reg28*reg138; reg122=reg116+reg122; reg55=reg77*reg55; reg126=reg28*reg126; reg77=reg77*reg144;
   reg74=reg104+reg74; reg144=reg54*reg144; reg30=reg117+reg30; reg138=reg129+reg138; reg9=reg51*reg141;
   reg141=reg121*reg141; reg111=reg114+reg111; reg126=reg125+reg126; reg1=reg51*reg1; reg120=reg120*reg140;
   reg119=reg119*reg139; reg148=reg147/reg148; reg140=reg85*reg140; reg139=reg86*reg139; reg77=reg74+reg77;
   reg115=reg85*reg115; reg128=reg86*reg128; reg145=reg133+reg145; reg55=reg122+reg55; reg144=reg134+reg144;
   reg24=reg143+reg24; reg77=reg124*reg77; reg1=reg126+reg1; reg119=reg120+reg119; reg128=reg115+reg128;
   reg108=reg68*reg108; reg127=reg127*reg148; reg24=reg131*reg24; reg144=reg55*reg144; reg141=reg111+reg141;
   reg139=reg140+reg139; reg145=reg30*reg145; reg148=reg68*reg148; reg9=reg138+reg9; reg127=reg119+reg127;
   reg148=reg139+reg148; reg24=reg145-reg24; reg108=reg128+reg108; reg77=reg144-reg77; reg141=reg1*reg141;
   reg9=reg132*reg9; reg1=0.024056261216234395431*reg77; reg10=0.024056261216234395431*reg24; reg28=0.035220810900864524453*reg24; reg30=0.035220810900864524453*reg77;
   reg31=0.024056261216234409915*reg77; reg9=reg141-reg9; reg37=0.04166666666666666908*reg24; reg38=0.13144585576580215187*reg24; reg148=reg135*reg148;
   reg51=0.13144585576580215187*reg77; reg54=0.04166666666666666908*reg77; reg127=reg108*reg127; reg55=0.024056261216234409915*reg24; reg61=0.13144585576580215187*reg9;
   reg68=reg38+reg51; reg38=reg30+reg38; reg30=reg28+reg30; reg69=0.035220810900864524453*reg9; reg74=0.024056261216234409915*reg9;
   reg1=reg1-reg37; reg75=0.04166666666666666908*reg9; reg55=reg55+reg54; reg148=reg127-reg148; reg80=0.024056261216234395431*reg9;
   reg31=reg37+reg31; reg51=reg28+reg51; reg54=reg10-reg54; reg55=reg55+reg75; reg10=0.024056261216234409915*reg148;
   reg28=0.035220810900864524453*reg148; reg38=reg61+reg38; reg75=reg54-reg75; reg37=0.024056261216234395431*reg148; reg68=reg68+reg69;
   reg54=0.13144585576580215187*reg148; reg61=reg30+reg61; reg74=reg1-reg74; reg31=reg80-reg31; reg1=0.04166666666666666908*reg148;
   reg51=reg69+reg51; Ne(0,21)+=reg54+reg51; Ne(1,22)+=reg54+reg51; Ne(2,23)+=reg54+reg51; Ne(0,18)+=reg28+reg68;
   Ne(1,19)+=reg28+reg68; Ne(2,20)+=reg28+reg68; Ne(0,15)+=reg38+reg28; Ne(1,16)+=reg38+reg28; Ne(2,17)+=reg38+reg28;
   Ne(0,12)+=reg61+reg54; Ne(1,13)+=reg61+reg54; Ne(2,14)+=reg61+reg54; Ne(0,9)+=reg74-reg1; Ne(1,10)+=reg74-reg1;
   Ne(2,11)+=reg74-reg1; Ne(0,6)+=reg75-reg10; Ne(1,7)+=reg75-reg10; Ne(2,8)+=reg75-reg10; Ne(0,3)+=reg31-reg1;
   Ne(1,4)+=reg31-reg1; Ne(2,5)+=reg31-reg1; Ne(0,0)+=reg37-reg55; Ne(1,1)+=reg37-reg55; Ne(2,2)+=reg37-reg55;

}

};
} // namespace LMT

