
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
   
   double reg0=0.21132486540518713447*e.pos(0)[0]; double reg1=0.21132486540518713447*e.pos(1)[0]; double reg2=0.21132486540518713447*e.pos(0)[1]; double reg3=0.21132486540518713447*e.pos(1)[1]; double reg4=0.78867513459481286553*e.pos(0)[0];
   double reg5=reg3-reg2; double reg6=0.78867513459481286553*e.pos(2)[1]; double reg7=0.21132486540518713447*e.pos(0)[2]; double reg8=0.21132486540518713447*e.pos(1)[2]; double reg9=0.78867513459481286553*e.pos(0)[1];
   double reg10=0.78867513459481286553*e.pos(2)[0]; double reg11=reg1-reg0; double reg12=0.78867513459481286553*e.pos(1)[1]; double reg13=0.78867513459481286553*e.pos(1)[0]; double reg14=0.21132486540518713447*e.pos(2)[1];
   double reg15=0.78867513459481286553*e.pos(2)[2]; double reg16=reg8-reg7; double reg17=0.21132486540518713447*e.pos(2)[0]; double reg18=0.78867513459481286553*e.pos(1)[2]; double reg19=0.78867513459481286553*e.pos(0)[2];
   double reg20=0.78867513459481286553*e.pos(3)[1]; reg5=reg5+reg6; double reg21=0.78867513459481286553*e.pos(3)[0]; reg11=reg11+reg10; double reg22=reg12-reg9;
   double reg23=reg13-reg4; double reg24=0.21132486540518713447*e.pos(2)[2]; double reg25=reg18-reg19; reg22=reg22+reg14; reg23=reg17+reg23;
   double reg26=0.21132486540518713447*e.pos(3)[1]; double reg27=0.21132486540518713447*e.pos(3)[0]; reg11=reg11-reg21; reg5=reg5-reg20; reg16=reg16+reg15;
   double reg28=0.78867513459481286553*e.pos(3)[2]; double reg29=pow(reg5,2); reg23=reg23-reg27; reg22=reg22-reg26; reg16=reg16-reg28;
   double reg30=pow(reg11,2); reg25=reg24+reg25; double reg31=0.21132486540518713447*e.pos(3)[2]; double reg32=pow(reg16,2); reg25=reg25-reg31;
   reg29=reg30+reg29; reg30=pow(reg22,2); double reg33=pow(reg23,2); reg9=reg3+reg9; reg3=pow(reg25,2);
   reg32=reg29+reg32; reg4=reg1+reg4; reg12=reg2+reg12; reg30=reg33+reg30; reg13=reg0+reg13;
   reg4=reg17-reg4; reg9=reg14-reg9; reg19=reg8+reg19; reg18=reg7+reg18; reg3=reg30+reg3;
   reg12=reg6-reg12; reg32=pow(reg32,0.5); reg13=reg10-reg13; reg9=reg20+reg9; reg19=reg24-reg19;
   reg4=reg21+reg4; reg3=pow(reg3,0.5); reg0=reg5/reg32; reg1=reg11/reg32; reg27=reg13+reg27;
   reg26=reg12+reg26; reg18=reg15-reg18; reg32=reg16/reg32; reg2=reg1*reg27; reg6=reg0*reg26;
   reg19=reg28+reg19; reg31=reg18+reg31; reg7=reg0*reg9; reg8=reg1*reg4; reg10=reg23/reg3;
   reg12=reg22/reg3; reg6=reg2+reg6; reg2=reg32*reg19; reg13=reg32*reg31; reg7=reg8+reg7;
   reg8=reg27*reg10; reg3=reg25/reg3; reg14=reg26*reg12; reg15=reg9*reg12; reg17=reg4*reg10;
   reg2=reg7+reg2; reg7=reg31*reg3; reg13=reg6+reg13; reg14=reg8+reg14; reg6=reg0*reg2;
   reg8=reg1*reg2; reg7=reg14+reg7; reg15=reg17+reg15; reg14=reg19*reg3; reg17=reg0*reg13;
   reg18=reg1*reg13; reg18=reg27-reg18; reg20=reg12*reg7; reg21=reg32*reg2; reg6=reg9-reg6;
   reg24=reg32*reg13; reg8=reg4-reg8; reg14=reg15+reg14; reg17=reg26-reg17; reg15=reg10*reg7;
   reg28=pow(reg18,2); reg29=reg12*reg14; reg30=reg10*reg14; reg33=reg3*reg7; reg20=reg26-reg20;
   reg15=reg27-reg15; reg21=reg19-reg21; reg24=reg31-reg24; double reg34=pow(reg6,2); double reg35=pow(reg17,2);
   double reg36=pow(reg8,2); double reg37=pow(reg20,2); double reg38=pow(reg15,2); reg33=reg31-reg33; double reg39=reg3*reg14;
   reg30=reg4-reg30; reg29=reg9-reg29; double reg40=pow(reg24,2); reg35=reg28+reg35; reg34=reg36+reg34;
   reg28=pow(reg21,2); reg36=pow(reg29,2); reg40=reg35+reg40; reg39=reg19-reg39; reg35=pow(reg33,2);
   reg28=reg34+reg28; reg37=reg38+reg37; reg34=pow(reg30,2); reg36=reg34+reg36; reg34=pow(reg39,2);
   reg28=pow(reg28,0.5); reg35=reg37+reg35; reg40=pow(reg40,0.5); reg6=reg6/reg28; reg8=reg8/reg28;
   reg17=reg17/reg40; reg35=pow(reg35,0.5); reg34=reg36+reg34; reg18=reg18/reg40; reg40=reg24/reg40;
   reg24=reg26*reg17; reg36=reg11*reg18; reg18=reg27*reg18; reg15=reg15/reg35; reg17=reg5*reg17;
   reg20=reg20/reg35; reg0=reg5*reg0; reg37=reg11*reg8; reg5=reg5*reg6; reg34=pow(reg34,0.5);
   reg28=reg21/reg28; reg6=reg9*reg6; reg8=reg4*reg8; reg1=reg11*reg1; reg12=reg22*reg12;
   reg11=reg22*reg20; reg21=reg23*reg15; reg5=reg37+reg5; reg37=reg16*reg28; reg35=reg33/reg35;
   reg20=reg26*reg20; reg15=reg27*reg15; reg29=reg29/reg34; reg30=reg30/reg34; reg0=reg1+reg0;
   reg28=reg19*reg28; reg1=reg16*reg40; reg24=reg18+reg24; reg6=reg8+reg6; reg40=reg31*reg40;
   reg32=reg16*reg32; reg17=reg36+reg17; reg10=reg23*reg10; reg4=reg4*reg30; reg9=reg9*reg29;
   reg34=reg39/reg34; reg30=reg23*reg30; reg32=reg0+reg32; reg29=reg22*reg29; reg28=reg6+reg28;
   reg0=reg25*reg35; reg11=reg21+reg11; reg37=reg5+reg37; reg12=reg10+reg12; reg35=reg31*reg35;
   reg20=reg15+reg20; reg40=reg24+reg40; reg1=reg17+reg1; reg3=reg25*reg3; reg9=reg4+reg9;
   reg28=reg32*reg28; reg19=reg19*reg34; reg37=reg2*reg37; reg1=reg13*reg1; reg3=reg12+reg3;
   reg40=reg32*reg40; reg35=reg20+reg35; reg29=reg30+reg29; reg34=reg25*reg34; reg0=reg11+reg0;
   reg19=reg9+reg19; reg34=reg29+reg34; reg37=reg28-reg37; reg1=reg40-reg1; reg35=reg3*reg35;
   reg0=reg7*reg0; reg2=0.1555021169820365473*reg37; reg4=0.1555021169820365473*reg1; reg5=0.011164549684630114537*reg37; reg6=0.04166666666666666908*reg1;
   reg34=reg14*reg34; reg0=reg35-reg0; reg19=reg3*reg19; reg3=0.011164549684630114537*reg1; reg7=0.04166666666666666908*reg37;
   reg8=0.04166666666666666908*reg0; reg34=reg19-reg34; reg9=0.011164549684630114537*reg0; reg2=reg6+reg2; reg5=reg6+reg5;
   reg6=0.1555021169820365473*reg0; reg3=reg3+reg7; reg4=reg7+reg4; reg9=reg2+reg9; reg2=0.011164549684630114537*reg34;
   reg4=reg8+reg4; reg7=0.04166666666666666908*reg34; reg6=reg5+reg6; reg5=0.1555021169820365473*reg34; reg8=reg3+reg8;
   Ne(0,3)+=reg6+reg7; Ne(1,4)+=reg6+reg7; Ne(2,5)+=reg6+reg7; Ne(0,6)+=reg4+reg2; Ne(1,7)+=reg4+reg2;
   Ne(2,8)+=reg4+reg2; Ne(0,0)+=reg8+reg5; Ne(1,1)+=reg8+reg5; Ne(2,2)+=reg8+reg5; Ne(0,9)+=reg7+reg9;
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
   
   double reg0=0.56758792732771999991*e.pos(0)[0]; double reg1=0.78379396366385999995*e.pos(1)[0]; double reg2=0.78379396366385999995*e.pos(0)[0]; double reg3=0.78379396366385999995*e.pos(1)[1]; double reg4=0.56758792732771999991*e.pos(0)[1];
   double reg5=0.56758792732771999991*e.pos(1)[0]; double reg6=0.56758792732771999991*e.pos(1)[1]; double reg7=0.78379396366385999995*e.pos(0)[1]; double reg8=reg3+reg4; double reg9=reg6+reg7;
   double reg10=1.3513818909915799999*e.pos(3)[1]; double reg11=reg5+reg2; double reg12=1.3513818909915799999*e.pos(3)[0]; double reg13=0.78379396366385999995*e.pos(1)[2]; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg1+reg0; double reg16=0.56758792732771999991*e.pos(1)[2]; double reg17=0.78379396366385999995*e.pos(0)[2]; double reg18=2.2673902919218320001*e.pos(0)[1]; double reg19=1.3513818909915799999*e.pos(3)[2];
   double reg20=reg13+reg14; double reg21=0.63369514596091600003*e.pos(1)[0]; double reg22=0.63369514596091600003*e.pos(1)[1]; double reg23=reg10-reg9; double reg24=2.2673902919218320001*e.pos(0)[0];
   double reg25=1.78379396366386*e.pos(4)[1]; reg8=reg8-reg10; double reg26=reg16+reg17; double reg27=reg12-reg11; double reg28=1.78379396366386*e.pos(4)[0];
   reg15=reg15-reg12; double reg29=2.2673902919218320001*e.pos(1)[0]; reg3=reg3-reg7; double reg30=0.63369514596091600003*e.pos(0)[1]; double reg31=2.2673902919218320001*e.pos(1)[1];
   double reg32=reg21+reg24; double reg33=0.63369514596091600003*e.pos(0)[0]; double reg34=2.9010854378827480001*e.pos(3)[0]; double reg35=2.9010854378827480001*e.pos(3)[1]; double reg36=2.2673902919218320001*e.pos(0)[2];
   double reg37=0.63369514596091600003*e.pos(1)[2]; double reg38=reg22+reg18; double reg39=0.43241207267228000009*e.pos(4)[0]; double reg40=0.43241207267228000009*e.pos(4)[1]; reg15=reg15+reg28;
   reg1=reg1-reg2; reg27=reg28+reg27; double reg41=1.78379396366386*e.pos(5)[0]; double reg42=1.78379396366386*e.pos(4)[2]; reg20=reg20-reg19;
   reg23=reg25+reg23; double reg43=1.78379396366386*e.pos(5)[1]; reg8=reg8+reg25; double reg44=reg19-reg26; double reg45=0.366304854039084*e.pos(4)[0];
   reg44=reg42+reg44; double reg46=0.43241207267228000009*e.pos(4)[2]; double reg47=reg35-reg38; double reg48=reg34-reg32; reg1=reg1+reg39;
   reg23=reg23-reg43; reg27=reg27-reg41; reg13=reg13-reg17; double reg49=0.43241207267228000009*e.pos(5)[1]; reg3=reg40+reg3;
   double reg50=0.43241207267228000009*e.pos(5)[0]; reg41=reg15-reg41; reg43=reg8-reg43; reg20=reg20+reg42; reg8=1.78379396366386*e.pos(5)[2];
   reg15=0.63369514596091600003*e.pos(0)[2]; double reg51=2.2673902919218320001*e.pos(1)[2]; reg31=reg31+reg30; reg29=reg29+reg33; double reg52=0.366304854039084*e.pos(4)[1];
   double reg53=2.9010854378827480001*e.pos(3)[2]; double reg54=reg37+reg36; double reg55=reg53-reg54; reg48=reg48+reg45; reg51=reg51+reg15;
   double reg56=reg30-reg22; reg50=reg1-reg50; reg1=pow(reg23,2); double reg57=pow(reg41,2); double reg58=reg33-reg21;
   double reg59=0.366304854039084*e.pos(5)[0]; double reg60=2.710505431213761085e-20*e.pos(3)[0]; reg20=reg20-reg8; double reg61=pow(reg27,2); double reg62=pow(reg43,2);
   double reg63=2.710505431213761085e-20*e.pos(3)[1]; double reg64=0.78379396366385999995*e.pos(2)[1]; reg29=reg29-reg34; reg47=reg47+reg52; double reg65=0.366304854039084*e.pos(4)[2];
   double reg66=0.78379396366385999995*e.pos(2)[0]; reg31=reg31-reg35; reg8=reg44-reg8; reg44=0.366304854039084*e.pos(5)[1]; double reg67=0.43241207267228000009*e.pos(5)[2];
   reg13=reg46+reg13; reg49=reg3-reg49; reg47=reg47-reg44; reg3=reg15-reg37; double reg68=2.710505431213761085e-20*e.pos(3)[2];
   reg56=reg56-reg63; reg48=reg48-reg59; reg58=reg58-reg60; double reg69=3.2673902919218320001*e.pos(4)[1]; double reg70=3.2673902919218320001*e.pos(4)[0];
   reg51=reg51-reg53; double reg71=pow(reg50,2); double reg72=pow(reg49,2); reg31=reg52+reg31; reg67=reg13-reg67;
   reg13=0.56758792732771999991*e.pos(2)[0]; double reg73=0.56758792732771999991*e.pos(2)[1]; reg29=reg45+reg29; reg55=reg55+reg65; double reg74=0.366304854039084*e.pos(5)[2];
   double reg75=pow(reg20,2); double reg76=pow(reg8,2); double reg77=reg0+reg66; double reg78=1.78379396366386*e.pos(3)[0]; double reg79=1.78379396366386*e.pos(3)[1];
   double reg80=reg4+reg64; double reg81=0.78379396366385999995*e.pos(2)[2]; reg1=reg61+reg1; reg62=reg57+reg62; reg58=reg70+reg58;
   reg57=3.2673902919218320001*e.pos(5)[0]; reg56=reg69+reg56; reg61=3.2673902919218320001*e.pos(5)[1]; double reg82=0.63369514596091600003*e.pos(2)[1]; reg3=reg3-reg68;
   double reg83=0.63369514596091600003*e.pos(2)[0]; reg75=reg62+reg75; reg62=3.2673902919218320001*e.pos(4)[2]; reg59=reg29-reg59; reg29=0.56758792732771999991*e.pos(2)[2];
   reg77=reg77-reg78; reg72=reg71+reg72; reg71=reg7+reg73; double reg84=pow(reg67,2); double reg85=reg2+reg13;
   reg51=reg65+reg51; reg44=reg31-reg44; reg31=1.78379396366386*e.pos(3)[2]; double reg86=reg14+reg81; double reg87=pow(reg47,2);
   reg7=reg64-reg7; reg55=reg55-reg74; reg64=pow(reg48,2); reg2=reg66-reg2; reg66=0.43241207267228000009*e.pos(3)[1];
   reg80=reg80-reg79; double reg88=0.43241207267228000009*e.pos(3)[0]; reg76=reg1+reg76; reg1=0.63369514596091600003*e.pos(2)[2]; reg74=reg51-reg74;
   reg78=reg78+reg85; reg76=pow(reg76,0.5); reg77=reg28+reg77; reg88=reg2-reg88; reg2=reg17+reg29;
   reg51=pow(reg59,2); reg84=reg72+reg84; reg87=reg64+reg87; reg79=reg79+reg71; reg86=reg86-reg31;
   reg64=pow(reg44,2); reg24=reg24+reg83; reg72=0.366304854039084*e.pos(3)[0]; double reg89=3.2673902919218320001*e.pos(5)[2]; reg3=reg62+reg3;
   double reg90=pow(reg55,2); double reg91=0.43241207267228000009*e.pos(3)[2]; reg61=reg56-reg61; reg17=reg81-reg17; reg80=reg25+reg80;
   reg18=reg18+reg82; reg57=reg58-reg57; reg56=1.3513818909915799999*e.pos(5)[1]; reg66=reg7-reg66; reg75=pow(reg75,0.5);
   reg7=0.366304854039084*e.pos(3)[1]; reg58=1.3513818909915799999*e.pos(5)[0]; reg28=reg28-reg78; reg86=reg42+reg86; reg81=1.3513818909915799999*e.pos(5)[2];
   reg25=reg25-reg79; reg31=reg31+reg2; reg77=reg77-reg58; reg80=reg80-reg56; double reg92=0.366304854039084*e.pos(3)[2];
   reg36=reg36+reg1; double reg93=reg18+reg7; double reg94=reg24+reg72; reg90=reg87+reg90; reg87=reg27/reg76;
   double reg95=reg23/reg76; double reg96=2.2673902919218320001*e.pos(2)[1]; reg39=reg88+reg39; reg40=reg66+reg40; reg91=reg17-reg91;
   reg17=2.2673902919218320001*e.pos(2)[0]; reg89=reg3-reg89; reg3=pow(reg61,2); reg66=pow(reg57,2); reg88=reg41/reg75;
   double reg97=3.2673902919218320001*e.pos(3)[1]; double reg98=reg30-reg82; double reg99=reg43/reg75; reg64=reg51+reg64; reg51=3.2673902919218320001*e.pos(3)[0];
   double reg100=reg33-reg83; double reg101=pow(reg74,2); reg84=pow(reg84,0.5); reg3=reg66+reg3; reg101=reg64+reg101;
   reg64=pow(reg89,2); reg66=2.9010854378827480001*e.pos(5)[0]; reg90=pow(reg90,0.5); reg17=reg33+reg17; reg46=reg91+reg46;
   reg28=reg58+reg28; reg91=reg95*reg40; double reg102=reg87*reg39; double reg103=reg99*reg80; double reg104=reg50/reg84;
   double reg105=reg49/reg84; reg96=reg30+reg96; reg76=reg8/reg76; double reg106=2.2673902919218320001*e.pos(2)[2]; reg75=reg20/reg75;
   reg86=reg86-reg81; reg51=reg100-reg51; reg42=reg42-reg31; reg25=reg56+reg25; reg97=reg98-reg97;
   reg98=reg88*reg77; reg100=reg36+reg92; double reg107=reg15-reg1; double reg108=3.2673902919218320001*e.pos(3)[2]; double reg109=2.9010854378827480001*e.pos(5)[1];
   double reg110=reg45-reg94; double reg111=reg52-reg93; double reg112=reg75*reg86; reg70=reg51+reg70; reg106=reg15+reg106;
   reg110=reg66+reg110; reg101=pow(reg101,0.5); reg7=reg96-reg7; reg108=reg107-reg108; reg64=reg3+reg64;
   reg69=reg97+reg69; reg3=2.9010854378827480001*e.pos(5)[2]; reg103=reg98+reg103; reg51=reg65-reg100; reg96=5.42101086242752217e-20*e.pos(5)[0];
   reg97=5.42101086242752217e-20*e.pos(5)[1]; reg98=reg47/reg90; reg72=reg17-reg72; reg17=reg48/reg90; reg111=reg111+reg109;
   reg91=reg102+reg91; reg102=reg76*reg46; reg84=reg67/reg84; reg107=reg104*reg28; double reg113=reg105*reg25;
   reg42=reg81+reg42; reg72=reg45+reg72; reg7=reg52+reg7; reg45=reg17*reg110; reg102=reg91+reg102;
   reg90=reg55/reg90; reg52=reg98*reg111; reg92=reg106-reg92; reg51=reg51+reg3; reg64=pow(reg64,0.5);
   reg91=reg84*reg42; reg70=reg70-reg96; reg112=reg103+reg112; reg69=reg69-reg97; reg62=reg108+reg62;
   reg103=5.42101086242752217e-20*e.pos(5)[2]; reg106=reg44/reg101; reg108=reg59/reg101; reg113=reg107+reg113; reg107=reg61/reg64;
   double reg114=reg90*reg51; double reg115=reg57/reg64; reg52=reg45+reg52; reg101=reg74/reg101; reg62=reg62-reg103;
   reg45=reg87*reg102; double reg116=reg106*reg69; double reg117=reg95*reg102; double reg118=reg108*reg70; reg91=reg113+reg91;
   reg113=reg88*reg112; double reg119=reg99*reg112; reg92=reg65+reg92; reg72=reg72-reg66; reg7=reg7-reg109;
   reg113=reg77-reg113; reg65=reg115*reg72; reg116=reg118+reg116; reg119=reg80-reg119; reg118=reg101*reg62;
   double reg120=reg75*reg112; double reg121=reg104*reg91; reg117=reg40-reg117; double reg122=reg76*reg102; reg45=reg39-reg45;
   reg64=reg89/reg64; reg114=reg52+reg114; reg92=reg92-reg3; reg52=reg107*reg7; double reg123=reg105*reg91;
   double reg124=pow(reg45,2); double reg125=reg84*reg91; reg122=reg46-reg122; reg121=reg28-reg121; double reg126=pow(reg117,2);
   reg120=reg86-reg120; reg118=reg116+reg118; reg116=pow(reg119,2); double reg127=reg98*reg114; double reg128=pow(reg113,2);
   reg52=reg65+reg52; reg65=reg17*reg114; reg123=reg25-reg123; double reg129=reg64*reg92; reg116=reg128+reg116;
   reg128=reg108*reg118; double reg130=pow(reg120,2); double reg131=reg106*reg118; reg125=reg42-reg125; double reg132=pow(reg122,2);
   reg124=reg126+reg124; reg126=pow(reg121,2); reg65=reg110-reg65; reg127=reg111-reg127; reg129=reg52+reg129;
   reg52=reg90*reg114; double reg133=pow(reg123,2); double reg134=pow(reg65,2); double reg135=pow(reg127,2); double reg136=pow(reg125,2);
   reg52=reg51-reg52; reg128=reg70-reg128; reg131=reg69-reg131; double reg137=reg101*reg118; reg130=reg116+reg130;
   reg133=reg126+reg133; reg116=reg115*reg129; reg126=reg107*reg129; reg132=reg124+reg132; reg116=reg72-reg116;
   reg124=reg64*reg129; reg136=reg133+reg136; reg135=reg134+reg135; reg132=pow(reg132,0.5); reg133=pow(reg52,2);
   reg134=pow(reg128,2); double reg138=pow(reg131,2); reg137=reg62-reg137; reg126=reg7-reg126; reg130=pow(reg130,0.5);
   double reg139=pow(reg116,2); double reg140=pow(reg126,2); reg119=reg119/reg130; reg136=pow(reg136,0.5); reg113=reg113/reg130;
   reg124=reg92-reg124; reg45=reg45/reg132; reg133=reg135+reg133; reg117=reg117/reg132; reg138=reg134+reg138;
   reg134=pow(reg137,2); reg39=reg39*reg45; reg45=reg27*reg45; reg132=reg122/reg132; reg40=reg40*reg117;
   reg117=reg23*reg117; reg140=reg139+reg140; reg122=pow(reg124,2); reg95=reg23*reg95; reg87=reg27*reg87;
   reg88=reg41*reg88; reg134=reg138+reg134; reg99=reg43*reg99; reg130=reg120/reg130; reg41=reg41*reg113;
   reg43=reg43*reg119; reg133=pow(reg133,0.5); reg121=reg121/reg136; reg123=reg123/reg136; reg113=reg77*reg113;
   reg119=reg80*reg119; reg99=reg88+reg99; reg75=reg20*reg75; reg134=pow(reg134,0.5); reg117=reg45+reg117;
   reg122=reg140+reg122; reg46=reg46*reg132; reg40=reg39+reg40; reg119=reg113+reg119; reg23=reg49*reg123;
   reg27=reg50*reg121; reg127=reg127/reg133; reg65=reg65/reg133; reg95=reg87+reg95; reg76=reg8*reg76;
   reg121=reg28*reg121; reg136=reg125/reg136; reg123=reg25*reg123; reg132=reg8*reg132; reg43=reg41+reg43;
   reg20=reg20*reg130; reg130=reg86*reg130; reg105=reg49*reg105; reg104=reg50*reg104; reg84=reg67*reg84;
   reg46=reg40+reg46; reg75=reg99+reg75; reg8=reg47*reg127; reg25=reg48*reg65; reg133=reg52/reg133;
   reg127=reg111*reg127; reg65=reg110*reg65; reg132=reg117+reg132; reg67=reg67*reg136; reg20=reg43+reg20;
   reg123=reg121+reg123; reg136=reg42*reg136; reg76=reg95+reg76; reg98=reg47*reg98; reg17=reg48*reg17;
   reg23=reg27+reg23; reg105=reg104+reg105; reg122=pow(reg122,0.5); reg130=reg119+reg130; reg128=reg128/reg134;
   reg131=reg131/reg134; reg126=reg126/reg122; reg51=reg51*reg133; reg127=reg65+reg127; reg27=reg44*reg131;
   reg28=reg59*reg128; reg132=reg102*reg132; reg134=reg137/reg134; reg131=reg69*reg131; reg20=reg112*reg20;
   reg136=reg123+reg136; reg90=reg55*reg90; reg98=reg17+reg98; reg128=reg70*reg128; reg67=reg23+reg67;
   reg106=reg44*reg106; reg116=reg116/reg122; reg46=reg76*reg46; reg84=reg105+reg84; reg108=reg59*reg108;
   reg130=reg75*reg130; reg133=reg55*reg133; reg8=reg25+reg8; reg67=reg91*reg67; reg20=reg130-reg20;
   reg72=reg72*reg116; reg90=reg98+reg90; reg136=reg84*reg136; reg106=reg108+reg106; reg101=reg74*reg101;
   reg116=reg57*reg116; reg7=reg7*reg126; reg51=reg127+reg51; reg74=reg74*reg134; reg133=reg8+reg133;
   reg27=reg28+reg27; reg115=reg57*reg115; reg107=reg61*reg107; reg132=reg46-reg132; reg122=reg124/reg122;
   reg134=reg62*reg134; reg126=reg61*reg126; reg131=reg128+reg131; reg8=0.088847818743090689935*reg20; reg126=reg116+reg126;
   reg7=reg72+reg7; reg17=reg89*reg122; reg23=0.005384432036113586778*reg20; reg25=0.088847818743090689935*reg132; reg122=reg92*reg122;
   reg28=0.021537728144454347112*reg132; reg39=0.021537728144454347112*reg20; reg40=0.009463616120767210603*reg132; reg64=reg89*reg64; reg107=reg115+reg107;
   reg41=0.005384432036113586778*reg132; reg42=0.009463616120767210603*reg20; reg133=reg114*reg133; reg67=reg136-reg67; reg51=reg90*reg51;
   reg74=reg27+reg74; reg134=reg131+reg134; reg101=reg106+reg101; reg25=reg39+reg25; reg40=reg23+reg40;
   reg23=reg41+reg23; reg27=0.009463616120767210603*reg67; reg43=0.021537728144454347112*reg67; reg8=reg28+reg8; reg41=reg42+reg41;
   reg42=0.005384432036113586778*reg67; reg28=reg39+reg28; reg39=0.088847818743090689935*reg67; reg64=reg107+reg64; reg17=reg126+reg17;
   reg133=reg51-reg133; reg74=reg118*reg74; reg134=reg101*reg134; reg122=reg7+reg122; reg7=0.016449618187943419918*reg133;
   reg39=reg28+reg39; reg74=reg134-reg74; reg41=reg41+reg42; reg122=reg64*reg122; reg27=reg23+reg27;
   reg8=reg8+reg43; reg23=0.028457289286966203713*reg133; reg28=0.0018441552587796664112*reg133; reg44=0.0041124045469858549794*reg133; reg40=reg42+reg40;
   reg17=reg129*reg17; reg25=reg43+reg25; reg25=reg7+reg25; reg28=reg8+reg28; reg8=0.0018441552587796664109*reg74;
   reg42=0.016449618187943419918*reg74; reg43=0.016449618187943419916*reg74; reg45=0.004112404546985854979*reg74; reg41=reg23-reg41; reg17=reg122-reg17;
   reg23=0.028457289286966203713*reg74; reg40=reg40+reg44; reg27=reg44+reg27; reg44=0.0041124045469858549794*reg74; reg7=reg39+reg7;
   reg39=0.028457289286966203713*reg17; reg44=reg27+reg44; reg40=reg23-reg40; reg8=reg25+reg8; reg23=0.0018441552587796664111*reg17;
   reg42=reg28+reg42; reg25=0.016449618187943419918*reg17; reg45=reg41-reg45; reg27=0.0041124045469858549794*reg17; reg43=reg7+reg43;
   Ne(0,15)+=reg25+reg8; Ne(1,16)+=reg25+reg8; Ne(2,17)+=reg25+reg8; Ne(0,6)+=reg39-reg44; Ne(1,7)+=reg39-reg44;
   Ne(2,8)+=reg39-reg44; Ne(0,0)+=reg45-reg27; Ne(1,1)+=reg45-reg27; Ne(2,2)+=reg45-reg27; Ne(0,12)+=reg42+reg25;
   Ne(1,13)+=reg42+reg25; Ne(2,14)+=reg42+reg25; Ne(0,9)+=reg43+reg23; Ne(1,10)+=reg43+reg23; Ne(2,11)+=reg43+reg23;
   Ne(0,3)+=reg40-reg27; Ne(1,4)+=reg40-reg27; Ne(2,5)+=reg40-reg27;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.12200846792814621817*e.pos(1)[1]; double reg1=0.36602540378443865451*e.pos(0)[1]; double reg2=0.36602540378443865451*e.pos(0)[0]; double reg3=0.12200846792814621817*e.pos(1)[0]; double reg4=0.12200846792814621817*e.pos(0)[1];
   double reg5=0.36602540378443865451*e.pos(1)[1]; double reg6=0.12200846792814621817*e.pos(0)[0]; double reg7=0.36602540378443865451*e.pos(1)[0]; double reg8=reg0+reg1; double reg9=1.3660254037844385386*e.pos(2)[1];
   double reg10=0.45534180126147951289*e.pos(0)[0]; double reg11=0.45534180126147951289*e.pos(0)[1]; double reg12=1.3660254037844385386*e.pos(1)[0]; double reg13=0.12200846792814621817*e.pos(0)[2]; double reg14=0.36602540378443865451*e.pos(1)[2];
   double reg15=reg3+reg2; double reg16=0.45534180126147951289*e.pos(2)[0]; double reg17=1.3660254037844385386*e.pos(2)[0]; double reg18=0.45534180126147951289*e.pos(2)[1]; double reg19=reg7+reg6;
   double reg20=reg5+reg4; double reg21=1.3660254037844385386*e.pos(1)[1]; double reg22=0.12200846792814621817*e.pos(1)[2]; double reg23=0.36602540378443865451*e.pos(0)[2]; double reg24=0.45534180126147951289*e.pos(1)[0];
   double reg25=0.45534180126147951289*e.pos(1)[1]; double reg26=1.3660254037844385386*e.pos(3)[0]; double reg27=1.3660254037844385386*e.pos(2)[2]; reg8=reg8+reg9; double reg28=reg22+reg23;
   double reg29=0.45534180126147951289*e.pos(3)[1]; double reg30=0.45534180126147951289*e.pos(2)[2]; double reg31=reg14+reg13; double reg32=reg16+reg19; double reg33=1.3660254037844385386*e.pos(3)[1];
   double reg34=reg20+reg18; double reg35=1.3660254037844385386*e.pos(1)[2]; double reg36=reg11+reg21; double reg37=1.3660254037844385386*e.pos(0)[1]; double reg38=reg10+reg12;
   reg15=reg17+reg15; double reg39=0.45534180126147951289*e.pos(3)[0]; double reg40=0.45534180126147951289*e.pos(0)[2]; double reg41=0.12200846792814621817*e.pos(2)[1]; double reg42=0.12200846792814621817*e.pos(2)[0];
   double reg43=1.3660254037844385386*e.pos(0)[0]; double reg44=0.45534180126147951289*e.pos(1)[2]; double reg45=reg24+reg43; double reg46=reg25+reg37; double reg47=reg32+reg26;
   double reg48=0.36602540378443865451*e.pos(2)[0]; double reg49=reg40+reg35; double reg50=0.36602540378443865451*e.pos(3)[1]; reg36=reg41+reg36; double reg51=0.36602540378443865451*e.pos(3)[0];
   reg38=reg42+reg38; double reg52=reg34+reg33; double reg53=0.36602540378443865451*e.pos(2)[1]; double reg54=1.3660254037844385386*e.pos(0)[2]; double reg55=0.12200846792814621817*e.pos(2)[2];
   double reg56=1.3660254037844385386*e.pos(3)[2]; double reg57=reg31+reg30; reg15=reg15+reg39; double reg58=0.48803387171258487271*e.pos(4)[0]; double reg59=0.45534180126147951289*e.pos(3)[2];
   reg28=reg28+reg27; reg8=reg8+reg29; double reg60=0.48803387171258487271*e.pos(4)[1]; double reg61=1.8213672050459180516*e.pos(4)[1]; reg36=reg36+reg50;
   double reg62=1.8213672050459180516*e.pos(4)[0]; reg38=reg38+reg51; reg8=reg8-reg60; double reg63=reg46+reg53; double reg64=0.66666666666666670528*e.pos(5)[1];
   reg15=reg15-reg58; double reg65=reg44+reg54; double reg66=0.36602540378443865451*e.pos(2)[2]; double reg67=0.66666666666666670528*e.pos(5)[0]; double reg68=reg60-reg52;
   double reg69=0.48803387171258487271*e.pos(4)[2]; reg28=reg28+reg59; double reg70=reg57+reg56; double reg71=reg58-reg47; double reg72=0.12200846792814621817*e.pos(3)[0];
   double reg73=0.36602540378443865451*e.pos(3)[2]; reg49=reg55+reg49; double reg74=0.12200846792814621817*e.pos(3)[1]; double reg75=reg48+reg45; reg49=reg49+reg73;
   double reg76=1.8213672050459180516*e.pos(4)[2]; double reg77=reg74+reg63; double reg78=reg65+reg66; reg15=reg15+reg67; double reg79=1.8213672050459180516*e.pos(6)[0];
   double reg80=0.66666666666666670528*e.pos(5)[2]; double reg81=reg72+reg75; reg28=reg28-reg69; reg71=reg67+reg71; double reg82=reg69-reg70;
   double reg83=1.8213672050459180516*e.pos(6)[1]; reg8=reg8+reg64; reg68=reg64+reg68; double reg84=0.12200846792814621817*e.pos(3)[2]; reg36=reg36-reg61;
   reg38=reg38-reg62; reg36=reg64+reg36; reg8=reg8-reg83; double reg85=0.66666666666666670528*e.pos(7)[1]; double reg86=0.48803387171258487271*e.pos(6)[1];
   reg82=reg80+reg82; reg71=reg79+reg71; double reg87=reg1+reg25; double reg88=1.8213672050459180516*e.pos(6)[2]; reg28=reg28+reg80;
   double reg89=0.66666666666666670528*e.pos(7)[0]; reg15=reg15-reg79; reg68=reg83+reg68; double reg90=reg7+reg10; double reg91=reg62-reg81;
   double reg92=reg5+reg11; reg38=reg67+reg38; double reg93=reg61-reg77; double reg94=0.48803387171258487271*e.pos(6)[0]; reg49=reg49-reg76;
   double reg95=reg84+reg78; double reg96=reg2+reg24; reg21=reg4+reg21; reg93=reg64+reg93; reg71=reg71-reg89;
   reg38=reg38-reg94; reg8=reg8-reg85; reg4=0.48803387171258487271*e.pos(6)[2]; reg49=reg80+reg49; reg91=reg67+reg91;
   reg36=reg36-reg86; reg68=reg68-reg85; reg15=reg15-reg89; reg9=reg9+reg87; reg64=reg14+reg40;
   reg67=0.66666666666666670528*e.pos(7)[2]; reg41=reg41+reg92; double reg97=reg76-reg95; double reg98=reg23+reg44; reg42=reg42+reg90;
   reg17=reg17+reg96; reg12=reg6+reg12; reg28=reg28-reg88; reg82=reg88+reg82; reg35=reg13+reg35;
   reg82=reg82-reg67; reg43=reg3+reg43; reg38=reg38-reg89; reg93=reg86+reg93; reg42=reg26+reg42;
   reg49=reg49-reg4; reg21=reg18+reg21; reg97=reg80+reg97; reg91=reg94+reg91; reg41=reg33+reg41;
   reg3=pow(reg15,2); reg55=reg55+reg64; reg36=reg36-reg85; reg12=reg16+reg12; reg37=reg0+reg37;
   reg27=reg27+reg98; reg0=0.66666666666666670528*e.pos(4)[1]; reg74=reg9+reg74; reg6=0.66666666666666670528*e.pos(4)[0]; reg17=reg72+reg17;
   reg9=pow(reg71,2); reg13=pow(reg8,2); reg28=reg28-reg67; reg26=pow(reg68,2); reg41=reg41-reg0;
   reg33=0.48803387171258487271*e.pos(5)[1]; reg54=reg22+reg54; reg55=reg56+reg55; reg97=reg4+reg97; reg85=reg93-reg85;
   reg22=0.66666666666666670528*e.pos(4)[2]; reg56=pow(reg38,2); reg84=reg27+reg84; reg27=1.8213672050459180516*e.pos(5)[1]; reg74=reg74-reg0;
   reg72=pow(reg36,2); reg80=pow(reg28,2); reg93=1.8213672050459180516*e.pos(5)[0]; reg17=reg17-reg6; reg49=reg49-reg67;
   reg35=reg30+reg35; reg89=reg91-reg89; reg91=reg51+reg12; double reg99=reg50+reg21; reg43=reg48+reg43;
   double reg100=pow(reg82,2); reg37=reg53+reg37; reg42=reg42-reg6; double reg101=0.48803387171258487271*e.pos(5)[0]; reg13=reg3+reg13;
   reg26=reg9+reg26; reg42=reg42-reg101; reg67=reg97-reg67; reg3=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg22;
   reg9=reg29+reg37; reg97=0.66666666666666670528*e.pos(6)[1]; reg74=reg74-reg27; double reg102=pow(reg85,2); reg72=reg56+reg72;
   reg56=0.66666666666666670528*e.pos(6)[0]; reg17=reg17-reg93; double reg103=reg39+reg43; double reg104=pow(reg49,2); double reg105=0.48803387171258487271*e.pos(5)[2];
   reg55=reg55-reg22; double reg106=pow(reg89,2); reg54=reg66+reg54; reg80=reg13+reg80; reg100=reg26+reg100;
   reg13=reg6+reg91; reg26=reg0+reg99; double reg107=reg73+reg35; reg41=reg41-reg33; reg102=reg106+reg102;
   reg106=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg3; reg80=pow(reg80,0.5); double reg108=reg22+reg107; reg6=reg6+reg103;
   double reg109=reg27-reg26; double reg110=0.48803387171258487271*e.pos(7)[1]; reg74=reg74+reg97; double reg111=reg93-reg13; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg17=reg17+reg56; reg100=pow(reg100,0.5); reg104=reg72+reg104; reg0=reg0+reg9; reg41=reg97+reg41;
   reg72=1.8213672050459180516*e.pos(7)[1]; double reg113=1.8213672050459180516*e.pos(7)[0]; double reg114=reg59+reg54; reg42=reg56+reg42; reg55=reg55-reg105;
   double reg115=pow(reg67,2); double reg116=reg3-reg108; double reg117=reg71/reg100; reg104=pow(reg104,0.5); reg55=reg106+reg55;
   double reg118=1.8213672050459180516*e.pos(7)[2]; double reg119=reg68/reg100; reg111=reg56+reg111; reg17=reg17-reg112; reg42=reg42-reg113;
   double reg120=reg101-reg6; double reg121=reg15/reg80; reg74=reg74-reg110; reg41=reg41-reg72; reg109=reg97+reg109;
   double reg122=reg8/reg80; reg115=reg102+reg115; reg102=reg33-reg0; reg22=reg22+reg114; reg84=reg84+reg106;
   double reg123=0.48803387171258487271*e.pos(7)[2]; reg116=reg106+reg116; reg80=reg28/reg80; double reg124=reg36/reg104; double reg125=reg38/reg104;
   reg102=reg97+reg102; reg55=reg55-reg118; reg97=reg119*reg41; reg84=reg84-reg123; reg111=reg112+reg111;
   double reg126=reg121*reg17; reg120=reg56+reg120; reg100=reg82/reg100; reg56=reg117*reg42; reg115=pow(reg115,0.5);
   double reg127=reg105-reg22; double reg128=reg122*reg74; reg109=reg110+reg109; reg102=reg72+reg102; reg128=reg126+reg128;
   reg104=reg49/reg104; reg127=reg106+reg127; reg116=reg123+reg116; reg106=reg89/reg115; reg97=reg56+reg97;
   reg56=reg124*reg109; reg126=reg80*reg84; reg120=reg113+reg120; double reg129=reg100*reg55; double reg130=reg85/reg115;
   double reg131=reg125*reg111; reg56=reg131+reg56; reg127=reg118+reg127; reg131=reg104*reg116; reg129=reg97+reg129;
   reg126=reg128+reg126; reg97=reg130*reg102; reg115=reg67/reg115; reg128=reg106*reg120; double reg132=reg121*reg126;
   double reg133=reg122*reg126; double reg134=reg115*reg127; double reg135=reg119*reg129; reg131=reg56+reg131; reg56=reg117*reg129;
   reg97=reg128+reg97; reg128=reg124*reg131; double reg136=reg125*reg131; double reg137=reg80*reg126; reg132=reg17-reg132;
   reg134=reg97+reg134; reg97=reg100*reg129; reg133=reg74-reg133; reg56=reg42-reg56; reg135=reg41-reg135;
   double reg138=pow(reg132,2); reg97=reg55-reg97; double reg139=reg104*reg131; double reg140=pow(reg135,2); reg128=reg109-reg128;
   reg136=reg111-reg136; double reg141=pow(reg56,2); double reg142=pow(reg133,2); double reg143=reg106*reg134; double reg144=reg130*reg134;
   reg137=reg84-reg137; reg144=reg102-reg144; reg143=reg120-reg143; double reg145=reg115*reg134; double reg146=pow(reg136,2);
   double reg147=pow(reg128,2); reg139=reg116-reg139; double reg148=pow(reg137,2); reg142=reg138+reg142; reg138=pow(reg97,2);
   reg140=reg141+reg140; reg147=reg146+reg147; reg141=pow(reg139,2); reg138=reg140+reg138; reg148=reg142+reg148;
   reg145=reg127-reg145; reg140=pow(reg143,2); reg142=pow(reg144,2); reg148=pow(reg148,0.5); reg142=reg140+reg142;
   reg138=pow(reg138,0.5); reg141=reg147+reg141; reg140=pow(reg145,2); reg132=reg132/reg148; reg133=reg133/reg148;
   reg56=reg56/reg138; reg141=pow(reg141,0.5); reg135=reg135/reg138; reg140=reg142+reg140; reg42=reg42*reg56;
   reg41=reg41*reg135; reg138=reg97/reg138; reg119=reg68*reg119; reg117=reg71*reg117; reg97=reg8*reg133;
   reg142=reg15*reg132; reg148=reg137/reg148; reg133=reg74*reg133; reg132=reg17*reg132; reg122=reg8*reg122;
   reg121=reg15*reg121; reg140=pow(reg140,0.5); reg136=reg136/reg141; reg128=reg128/reg141; reg135=reg68*reg135;
   reg56=reg71*reg56; reg84=reg84*reg148; reg133=reg132+reg133; reg97=reg142+reg97; reg148=reg28*reg148;
   reg8=reg36*reg128; reg15=reg38*reg136; reg80=reg28*reg80; reg122=reg121+reg122; reg124=reg36*reg124;
   reg119=reg117+reg119; reg141=reg139/reg141; reg100=reg82*reg100; reg128=reg109*reg128; reg55=reg55*reg138;
   reg41=reg42+reg41; reg135=reg56+reg135; reg138=reg82*reg138; reg144=reg144/reg140; reg125=reg38*reg125;
   reg143=reg143/reg140; reg136=reg111*reg136; reg120=reg120*reg143; reg116=reg116*reg141; reg102=reg102*reg144;
   reg140=reg145/reg140; reg80=reg122+reg80; reg143=reg89*reg143; reg128=reg136+reg128; reg8=reg15+reg8;
   reg141=reg49*reg141; reg144=reg85*reg144; reg55=reg41+reg55; reg138=reg135+reg138; reg106=reg89*reg106;
   reg100=reg119+reg100; reg130=reg85*reg130; reg124=reg125+reg124; reg148=reg97+reg148; reg104=reg49*reg104;
   reg84=reg133+reg84; reg102=reg120+reg102; reg138=reg129*reg138; reg84=reg80*reg84; reg127=reg127*reg140;
   reg55=reg100*reg55; reg130=reg106+reg130; reg115=reg67*reg115; reg144=reg143+reg144; reg140=reg67*reg140;
   reg116=reg128+reg116; reg148=reg126*reg148; reg141=reg8+reg141; reg104=reg124+reg104; reg140=reg144+reg140;
   reg127=reg102+reg127; reg138=reg55-reg138; reg148=reg84-reg148; reg115=reg130+reg115; reg141=reg131*reg141;
   reg116=reg104*reg116; reg8=0.024056261216234395431*reg148; reg15=0.13144585576580215187*reg138; reg17=0.024056261216234395431*reg138; reg28=0.024056261216234409915*reg138;
   reg36=0.13144585576580215187*reg148; reg38=0.04166666666666666908*reg148; reg41=0.035220810900864524453*reg148; reg42=0.024056261216234409915*reg148; reg49=0.035220810900864524453*reg138;
   reg140=reg134*reg140; reg127=reg115*reg127; reg141=reg116-reg141; reg55=0.04166666666666666908*reg138; reg56=reg41+reg49;
   reg67=0.024056261216234409915*reg141; reg17=reg17-reg38; reg68=0.13144585576580215187*reg141; reg49=reg49+reg36; reg8=reg8-reg55;
   reg36=reg36+reg15; reg71=0.035220810900864524453*reg141; reg28=reg38+reg28; reg38=0.024056261216234395431*reg141; reg15=reg41+reg15;
   reg140=reg127-reg140; reg41=0.04166666666666666908*reg141; reg55=reg42+reg55; reg42=0.035220810900864524453*reg140; reg49=reg68+reg49;
   reg74=0.13144585576580215187*reg140; reg36=reg36+reg71; reg15=reg71+reg15; reg68=reg56+reg68; reg67=reg17-reg67;
   reg17=0.024056261216234409915*reg140; reg8=reg8-reg41; reg56=0.04166666666666666908*reg140; reg28=reg38-reg28; reg38=0.024056261216234395431*reg140;
   reg41=reg55+reg41; Ne(0,9)+=reg67-reg56; Ne(1,10)+=reg67-reg56; Ne(2,11)+=reg67-reg56; Ne(0,12)+=reg68+reg74;
   Ne(1,13)+=reg68+reg74; Ne(2,14)+=reg68+reg74; Ne(0,6)+=reg8-reg17; Ne(1,7)+=reg8-reg17; Ne(2,8)+=reg8-reg17;
   Ne(0,15)+=reg49+reg42; Ne(1,16)+=reg49+reg42; Ne(2,17)+=reg49+reg42; Ne(0,3)+=reg28-reg56; Ne(1,4)+=reg28-reg56;
   Ne(2,5)+=reg28-reg56; Ne(0,0)+=reg38-reg41; Ne(1,1)+=reg38-reg41; Ne(2,2)+=reg38-reg41; Ne(0,18)+=reg42+reg36;
   Ne(1,19)+=reg42+reg36; Ne(2,20)+=reg42+reg36; Ne(0,21)+=reg74+reg15; Ne(1,22)+=reg74+reg15; Ne(2,23)+=reg74+reg15;

}

};
} // namespace LMT

