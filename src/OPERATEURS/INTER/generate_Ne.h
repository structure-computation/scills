
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
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.78379396366385999995*e.pos(1)[1]; double reg3=0.78379396366385999995*e.pos(0)[1]; double reg4=0.56758792732771999991*e.pos(1)[1];
   double reg5=0.56758792732771999991*e.pos(0)[1]; double reg6=0.78379396366385999995*e.pos(0)[0]; double reg7=0.56758792732771999991*e.pos(1)[0]; double reg8=reg4+reg3; double reg9=1.3513818909915799999*e.pos(3)[0];
   double reg10=0.78379396366385999995*e.pos(1)[2]; double reg11=reg0+reg1; double reg12=reg2+reg5; double reg13=1.3513818909915799999*e.pos(3)[1]; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=0.56758792732771999991*e.pos(1)[2]; double reg16=0.78379396366385999995*e.pos(0)[2]; double reg17=reg7+reg6; double reg18=reg10+reg14; double reg19=1.3513818909915799999*e.pos(3)[2];
   double reg20=reg13-reg8; double reg21=2.2673902919218320001*e.pos(0)[0]; double reg22=1.78379396366386*e.pos(4)[1]; double reg23=reg9-reg17; double reg24=0.63369514596091600003*e.pos(1)[0];
   double reg25=2.2673902919218320001*e.pos(0)[1]; double reg26=0.63369514596091600003*e.pos(1)[1]; double reg27=reg15+reg16; double reg28=1.78379396366386*e.pos(4)[0]; reg11=reg11-reg9;
   reg12=reg12-reg13; double reg29=2.9010854378827480001*e.pos(3)[0]; double reg30=2.2673902919218320001*e.pos(1)[1]; double reg31=reg24+reg21; double reg32=0.63369514596091600003*e.pos(0)[0];
   double reg33=0.63369514596091600003*e.pos(0)[1]; double reg34=2.2673902919218320001*e.pos(1)[0]; double reg35=2.2673902919218320001*e.pos(0)[2]; double reg36=0.63369514596091600003*e.pos(1)[2]; double reg37=reg26+reg25;
   reg2=reg2-reg3; double reg38=2.9010854378827480001*e.pos(3)[1]; reg11=reg28+reg11; double reg39=1.78379396366386*e.pos(5)[0]; reg0=reg0-reg6;
   reg12=reg12+reg22; double reg40=1.78379396366386*e.pos(5)[1]; reg20=reg22+reg20; reg18=reg18-reg19; double reg41=1.78379396366386*e.pos(4)[2];
   double reg42=reg19-reg27; reg23=reg28+reg23; double reg43=0.43241207267228000009*e.pos(4)[1]; double reg44=0.43241207267228000009*e.pos(4)[0]; reg42=reg41+reg42;
   double reg45=0.43241207267228000009*e.pos(4)[2]; double reg46=reg29-reg31; double reg47=reg38-reg37; reg20=reg20-reg40; reg23=reg23-reg39;
   double reg48=0.366304854039084*e.pos(4)[0]; reg10=reg10-reg16; double reg49=0.43241207267228000009*e.pos(5)[1]; reg2=reg43+reg2; reg0=reg44+reg0;
   double reg50=0.43241207267228000009*e.pos(5)[0]; reg39=reg11-reg39; reg40=reg12-reg40; reg18=reg18+reg41; reg11=1.78379396366386*e.pos(5)[2];
   reg12=0.63369514596091600003*e.pos(0)[2]; double reg51=2.2673902919218320001*e.pos(1)[2]; reg30=reg30+reg33; reg34=reg34+reg32; double reg52=0.366304854039084*e.pos(4)[1];
   double reg53=reg36+reg35; double reg54=2.9010854378827480001*e.pos(3)[2]; reg50=reg0-reg50; reg51=reg51+reg12; reg0=reg32-reg24;
   double reg55=0.366304854039084*e.pos(5)[0]; double reg56=2.710505431213761085e-20*e.pos(3)[0]; reg18=reg18-reg11; double reg57=pow(reg23,2); double reg58=reg33-reg26;
   double reg59=pow(reg40,2); double reg60=pow(reg39,2); double reg61=pow(reg20,2); double reg62=2.710505431213761085e-20*e.pos(3)[1]; double reg63=0.366304854039084*e.pos(5)[1];
   reg11=reg42-reg11; reg47=reg47+reg52; reg42=0.78379396366385999995*e.pos(2)[0]; double reg64=0.366304854039084*e.pos(4)[2]; double reg65=0.43241207267228000009*e.pos(5)[2];
   reg10=reg45+reg10; reg30=reg30-reg38; reg34=reg34-reg29; reg49=reg2-reg49; reg2=0.78379396366385999995*e.pos(2)[1];
   double reg66=reg54-reg53; reg46=reg48+reg46; reg66=reg66+reg64; double reg67=reg12-reg36; double reg68=0.366304854039084*e.pos(5)[2];
   double reg69=2.710505431213761085e-20*e.pos(3)[2]; reg58=reg58-reg62; reg34=reg48+reg34; reg47=reg47-reg63; double reg70=0.56758792732771999991*e.pos(2)[1];
   double reg71=0.56758792732771999991*e.pos(2)[0]; reg65=reg10-reg65; reg46=reg46-reg55; reg0=reg0-reg56; reg10=pow(reg49,2);
   reg30=reg52+reg30; double reg72=3.2673902919218320001*e.pos(4)[1]; double reg73=3.2673902919218320001*e.pos(4)[0]; reg51=reg51-reg54; double reg74=pow(reg50,2);
   double reg75=pow(reg18,2); double reg76=reg1+reg42; double reg77=pow(reg11,2); double reg78=1.78379396366386*e.pos(3)[0]; double reg79=1.78379396366386*e.pos(3)[1];
   double reg80=reg5+reg2; double reg81=0.78379396366385999995*e.pos(2)[2]; reg61=reg57+reg61; reg59=reg60+reg59; reg0=reg73+reg0;
   reg57=3.2673902919218320001*e.pos(5)[0]; reg58=reg72+reg58; reg60=3.2673902919218320001*e.pos(5)[1]; double reg82=0.63369514596091600003*e.pos(2)[1]; reg67=reg67-reg69;
   double reg83=0.63369514596091600003*e.pos(2)[0]; reg75=reg59+reg75; reg59=3.2673902919218320001*e.pos(4)[2]; double reg84=0.56758792732771999991*e.pos(2)[2]; reg55=reg34-reg55;
   reg34=reg3+reg70; reg10=reg74+reg10; reg76=reg76-reg78; reg74=pow(reg65,2); double reg85=reg6+reg71;
   reg51=reg64+reg51; reg63=reg30-reg63; reg30=reg14+reg81; double reg86=pow(reg47,2); reg3=reg2-reg3;
   reg66=reg66-reg68; reg2=pow(reg46,2); double reg87=1.78379396366386*e.pos(3)[2]; double reg88=0.43241207267228000009*e.pos(3)[1]; reg6=reg42-reg6;
   reg80=reg80-reg79; reg42=0.43241207267228000009*e.pos(3)[0]; reg77=reg61+reg77; reg61=0.63369514596091600003*e.pos(2)[2]; reg30=reg30-reg87;
   reg42=reg6-reg42; reg78=reg78+reg85; reg76=reg28+reg76; reg74=reg10+reg74; reg6=reg16+reg84;
   reg10=pow(reg55,2); reg79=reg79+reg34; reg77=pow(reg77,0.5); reg86=reg2+reg86; reg68=reg51-reg68;
   reg2=pow(reg63,2); reg21=reg21+reg83; reg51=0.366304854039084*e.pos(3)[0]; double reg89=3.2673902919218320001*e.pos(5)[2]; reg67=reg59+reg67;
   double reg90=pow(reg66,2); double reg91=0.43241207267228000009*e.pos(3)[2]; reg60=reg58-reg60; reg16=reg81-reg16; reg80=reg22+reg80;
   reg25=reg25+reg82; reg57=reg0-reg57; reg0=1.3513818909915799999*e.pos(5)[1]; reg88=reg3-reg88; reg75=pow(reg75,0.5);
   reg3=0.366304854039084*e.pos(3)[1]; reg58=1.3513818909915799999*e.pos(5)[0]; reg30=reg41+reg30; reg81=1.3513818909915799999*e.pos(5)[2]; reg28=reg28-reg78;
   reg22=reg22-reg79; reg87=reg87+reg6; reg76=reg76-reg58; reg80=reg80-reg0; double reg92=0.366304854039084*e.pos(3)[2];
   reg35=reg35+reg61; double reg93=reg25+reg3; double reg94=reg21+reg51; reg90=reg86+reg90; reg86=reg23/reg77;
   double reg95=reg20/reg77; double reg96=2.2673902919218320001*e.pos(2)[1]; reg42=reg44+reg42; reg43=reg88+reg43; reg91=reg16-reg91;
   reg16=2.2673902919218320001*e.pos(2)[0]; reg89=reg67-reg89; reg44=pow(reg60,2); reg67=pow(reg57,2); reg88=reg39/reg75;
   double reg97=3.2673902919218320001*e.pos(3)[1]; double reg98=reg33-reg82; double reg99=3.2673902919218320001*e.pos(3)[0]; reg2=reg10+reg2; reg74=pow(reg74,0.5);
   reg10=pow(reg68,2); double reg100=reg40/reg75; double reg101=reg32-reg83; reg44=reg67+reg44; reg10=reg2+reg10;
   reg2=pow(reg89,2); reg67=2.9010854378827480001*e.pos(5)[0]; reg90=pow(reg90,0.5); reg28=reg58+reg28; reg16=reg32+reg16;
   reg45=reg91+reg45; reg91=reg95*reg43; reg99=reg101-reg99; reg101=reg86*reg42; double reg102=reg100*reg80;
   double reg103=reg50/reg74; double reg104=reg49/reg74; reg96=reg33+reg96; reg77=reg11/reg77; double reg105=2.2673902919218320001*e.pos(2)[2];
   reg75=reg18/reg75; reg30=reg30-reg81; reg41=reg41-reg87; reg22=reg0+reg22; reg97=reg98-reg97;
   reg98=reg88*reg76; double reg106=reg35+reg92; double reg107=reg12-reg61; double reg108=3.2673902919218320001*e.pos(3)[2]; double reg109=2.9010854378827480001*e.pos(5)[1];
   double reg110=reg52-reg93; double reg111=reg48-reg94; double reg112=reg75*reg30; reg111=reg67+reg111; reg105=reg12+reg105;
   reg10=pow(reg10,0.5); reg2=reg44+reg2; reg3=reg96-reg3; reg72=reg97+reg72; reg44=2.9010854378827480001*e.pos(5)[2];
   reg96=reg64-reg106; reg102=reg98+reg102; reg97=5.42101086242752217e-20*e.pos(5)[1]; reg98=reg47/reg90; double reg113=reg46/reg90;
   reg110=reg110+reg109; reg51=reg16-reg51; reg99=reg73+reg99; reg16=5.42101086242752217e-20*e.pos(5)[0]; reg108=reg107-reg108;
   reg91=reg101+reg91; reg73=reg77*reg45; reg74=reg65/reg74; reg101=reg103*reg28; reg107=reg104*reg22;
   reg41=reg81+reg41; reg51=reg48+reg51; reg3=reg52+reg3; reg48=reg113*reg111; reg73=reg91+reg73;
   reg90=reg66/reg90; reg52=reg98*reg110; reg92=reg105-reg92; reg96=reg96+reg44; reg2=pow(reg2,0.5);
   reg91=reg74*reg41; reg99=reg99-reg16; reg112=reg102+reg112; reg72=reg72-reg97; reg59=reg108+reg59;
   reg102=5.42101086242752217e-20*e.pos(5)[2]; reg105=reg63/reg10; reg108=reg55/reg10; reg107=reg101+reg107; reg101=reg60/reg2;
   double reg114=reg90*reg96; double reg115=reg57/reg2; reg52=reg48+reg52; reg10=reg68/reg10; reg59=reg59-reg102;
   reg48=reg86*reg73; double reg116=reg105*reg72; double reg117=reg95*reg73; double reg118=reg108*reg99; reg91=reg107+reg91;
   reg107=reg88*reg112; double reg119=reg100*reg112; reg92=reg64+reg92; reg51=reg51-reg67; reg3=reg3-reg109;
   reg107=reg76-reg107; reg64=reg115*reg51; reg116=reg118+reg116; reg119=reg80-reg119; reg118=reg10*reg59;
   double reg120=reg75*reg112; reg117=reg43-reg117; double reg121=reg103*reg91; double reg122=reg77*reg73; reg48=reg42-reg48;
   reg2=reg89/reg2; reg114=reg52+reg114; reg92=reg92-reg44; reg52=reg101*reg3; double reg123=reg104*reg91;
   reg52=reg64+reg52; reg64=pow(reg107,2); double reg124=pow(reg119,2); reg118=reg116+reg118; reg120=reg30-reg120;
   reg116=pow(reg117,2); double reg125=reg74*reg91; reg122=reg45-reg122; double reg126=pow(reg48,2); reg123=reg22-reg123;
   reg121=reg28-reg121; double reg127=reg2*reg92; double reg128=reg113*reg114; double reg129=reg98*reg114; reg124=reg64+reg124;
   reg125=reg41-reg125; reg128=reg111-reg128; reg64=reg108*reg118; reg129=reg110-reg129; double reg130=pow(reg122,2);
   reg126=reg116+reg126; reg116=pow(reg120,2); double reg131=pow(reg123,2); double reg132=reg105*reg118; reg127=reg52+reg127;
   reg52=pow(reg121,2); double reg133=reg90*reg114; reg116=reg124+reg116; reg64=reg99-reg64; reg124=pow(reg128,2);
   reg132=reg72-reg132; double reg134=pow(reg125,2); double reg135=reg10*reg118; reg133=reg96-reg133; double reg136=reg115*reg127;
   reg130=reg126+reg130; reg126=reg101*reg127; double reg137=pow(reg129,2); reg131=reg52+reg131; reg52=pow(reg132,2);
   reg135=reg59-reg135; reg136=reg51-reg136; double reg138=pow(reg64,2); double reg139=pow(reg133,2); reg137=reg124+reg137;
   reg124=reg2*reg127; reg116=pow(reg116,0.5); reg130=pow(reg130,0.5); reg134=reg131+reg134; reg126=reg3-reg126;
   reg134=pow(reg134,0.5); reg139=reg137+reg139; reg124=reg92-reg124; reg131=pow(reg136,2); reg48=reg48/reg130;
   reg117=reg117/reg130; reg119=reg119/reg116; reg137=pow(reg126,2); reg107=reg107/reg116; double reg140=pow(reg135,2);
   reg52=reg138+reg52; reg43=reg43*reg117; reg137=reg131+reg137; reg130=reg122/reg130; reg122=reg23*reg48;
   reg86=reg23*reg86; reg117=reg20*reg117; reg48=reg42*reg48; reg140=reg52+reg140; reg95=reg20*reg95;
   reg88=reg39*reg88; reg20=pow(reg124,2); reg100=reg40*reg100; reg116=reg120/reg116; reg39=reg39*reg107;
   reg40=reg40*reg119; reg119=reg80*reg119; reg139=pow(reg139,0.5); reg121=reg121/reg134; reg107=reg76*reg107;
   reg123=reg123/reg134; reg100=reg88+reg100; reg23=reg49*reg123; reg104=reg49*reg104; reg75=reg18*reg75;
   reg140=pow(reg140,0.5); reg117=reg122+reg117; reg20=reg137+reg20; reg45=reg45*reg130; reg43=reg48+reg43;
   reg119=reg107+reg119; reg42=reg50*reg121; reg95=reg86+reg95; reg77=reg11*reg77; reg129=reg129/reg139;
   reg128=reg128/reg139; reg134=reg125/reg134; reg123=reg22*reg123; reg121=reg28*reg121; reg130=reg11*reg130;
   reg40=reg39+reg40; reg18=reg18*reg116; reg103=reg50*reg103; reg116=reg30*reg116; reg75=reg100+reg75;
   reg139=reg133/reg139; reg45=reg43+reg45; reg110=reg110*reg129; reg74=reg65*reg74; reg104=reg103+reg104;
   reg130=reg117+reg130; reg111=reg111*reg128; reg20=pow(reg20,0.5); reg116=reg119+reg116; reg18=reg40+reg18;
   reg123=reg121+reg123; reg77=reg95+reg77; reg128=reg46*reg128; reg41=reg41*reg134; reg129=reg47*reg129;
   reg98=reg47*reg98; reg113=reg46*reg113; reg23=reg42+reg23; reg132=reg132/reg140; reg134=reg65*reg134;
   reg64=reg64/reg140; reg96=reg96*reg139; reg129=reg128+reg129; reg139=reg66*reg139; reg110=reg111+reg110;
   reg130=reg73*reg130; reg18=reg112*reg18; reg41=reg123+reg41; reg90=reg66*reg90; reg98=reg113+reg98;
   reg134=reg23+reg134; reg11=reg63*reg132; reg126=reg126/reg20; reg116=reg75*reg116; reg99=reg99*reg64;
   reg132=reg72*reg132; reg64=reg55*reg64; reg136=reg136/reg20; reg74=reg104+reg74; reg140=reg135/reg140;
   reg108=reg55*reg108; reg105=reg63*reg105; reg45=reg77*reg45; reg11=reg64+reg11; reg41=reg74*reg41;
   reg90=reg98+reg90; reg20=reg124/reg20; reg59=reg59*reg140; reg132=reg99+reg132; reg134=reg91*reg134;
   reg18=reg116-reg18; reg139=reg129+reg139; reg22=reg60*reg126; reg23=reg57*reg136; reg96=reg110+reg96;
   reg10=reg68*reg10; reg105=reg108+reg105; reg115=reg57*reg115; reg101=reg60*reg101; reg130=reg45-reg130;
   reg136=reg51*reg136; reg126=reg3*reg126; reg140=reg68*reg140; reg126=reg136+reg126; reg3=0.088847818743090689935*reg18;
   reg28=0.009463616120767210603*reg130; reg92=reg92*reg20; reg30=0.088847818743090689935*reg130; reg22=reg23+reg22; reg23=0.021537728144454347112*reg130;
   reg39=0.021537728144454347112*reg18; reg40=0.005384432036113586778*reg18; reg20=reg89*reg20; reg2=reg89*reg2; reg101=reg115+reg101;
   reg42=0.005384432036113586778*reg130; reg43=0.009463616120767210603*reg18; reg139=reg114*reg139; reg134=reg41-reg134; reg96=reg90*reg96;
   reg140=reg11+reg140; reg59=reg132+reg59; reg10=reg105+reg10; reg30=reg39+reg30; reg28=reg40+reg28;
   reg40=reg42+reg40; reg11=0.009463616120767210603*reg134; reg41=0.021537728144454347112*reg134; reg3=reg23+reg3; reg42=reg43+reg42;
   reg43=0.005384432036113586778*reg134; reg23=reg39+reg23; reg39=0.088847818743090689935*reg134; reg2=reg101+reg2; reg20=reg22+reg20;
   reg139=reg96-reg139; reg140=reg118*reg140; reg59=reg10*reg59; reg92=reg126+reg92; reg10=0.016449618187943419918*reg139;
   reg39=reg23+reg39; reg140=reg59-reg140; reg42=reg42+reg43; reg92=reg2*reg92; reg11=reg40+reg11;
   reg3=reg3+reg41; reg2=0.028457289286966203713*reg139; reg22=0.0018441552587796664112*reg139; reg23=0.0041124045469858549794*reg139; reg28=reg43+reg28;
   reg20=reg127*reg20; reg30=reg41+reg30; reg30=reg10+reg30; reg22=reg3+reg22; reg3=0.0018441552587796664109*reg140;
   reg40=0.016449618187943419918*reg140; reg41=0.016449618187943419916*reg140; reg43=0.004112404546985854979*reg140; reg42=reg2-reg42; reg20=reg92-reg20;
   reg2=0.028457289286966203713*reg140; reg28=reg28+reg23; reg11=reg23+reg11; reg23=0.0041124045469858549794*reg140; reg10=reg39+reg10;
   reg39=0.028457289286966203713*reg20; reg23=reg11+reg23; reg28=reg2-reg28; reg3=reg30+reg3; reg2=0.0018441552587796664111*reg20;
   reg40=reg22+reg40; reg11=0.016449618187943419918*reg20; reg43=reg42-reg43; reg22=0.0041124045469858549794*reg20; reg41=reg10+reg41;
   Ne(0,15)+=reg11+reg3; Ne(1,16)+=reg11+reg3; Ne(2,17)+=reg11+reg3; Ne(0,6)+=reg39-reg23; Ne(1,7)+=reg39-reg23;
   Ne(2,8)+=reg39-reg23; Ne(0,0)+=reg43-reg22; Ne(1,1)+=reg43-reg22; Ne(2,2)+=reg43-reg22; Ne(0,12)+=reg40+reg11;
   Ne(1,13)+=reg40+reg11; Ne(2,14)+=reg40+reg11; Ne(0,9)+=reg41+reg2; Ne(1,10)+=reg41+reg2; Ne(2,11)+=reg41+reg2;
   Ne(0,3)+=reg28-reg22; Ne(1,4)+=reg28-reg22; Ne(2,5)+=reg28-reg22;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.12200846792814621817*e.pos(1)[0]; double reg1=0.36602540378443865451*e.pos(0)[0]; double reg2=0.12200846792814621817*e.pos(1)[1]; double reg3=0.36602540378443865451*e.pos(0)[1]; double reg4=0.12200846792814621817*e.pos(0)[0];
   double reg5=0.36602540378443865451*e.pos(1)[0]; double reg6=0.12200846792814621817*e.pos(0)[1]; double reg7=0.36602540378443865451*e.pos(1)[1]; double reg8=0.45534180126147951289*e.pos(0)[0]; double reg9=0.36602540378443865451*e.pos(0)[2];
   double reg10=0.45534180126147951289*e.pos(0)[1]; double reg11=0.12200846792814621817*e.pos(1)[2]; double reg12=reg2+reg3; double reg13=1.3660254037844385386*e.pos(2)[1]; double reg14=reg0+reg1;
   double reg15=1.3660254037844385386*e.pos(2)[0]; double reg16=0.12200846792814621817*e.pos(0)[2]; double reg17=1.3660254037844385386*e.pos(1)[0]; double reg18=0.36602540378443865451*e.pos(1)[2]; double reg19=reg5+reg4;
   double reg20=1.3660254037844385386*e.pos(1)[1]; double reg21=0.45534180126147951289*e.pos(2)[1]; double reg22=reg7+reg6; double reg23=0.45534180126147951289*e.pos(2)[0]; double reg24=reg18+reg16;
   double reg25=0.45534180126147951289*e.pos(3)[1]; reg12=reg12+reg13; double reg26=1.3660254037844385386*e.pos(3)[1]; double reg27=reg22+reg21; double reg28=0.45534180126147951289*e.pos(2)[2];
   double reg29=1.3660254037844385386*e.pos(0)[0]; double reg30=0.45534180126147951289*e.pos(1)[0]; double reg31=1.3660254037844385386*e.pos(1)[2]; double reg32=1.3660254037844385386*e.pos(0)[1]; double reg33=reg10+reg20;
   double reg34=reg23+reg19; double reg35=0.45534180126147951289*e.pos(1)[1]; double reg36=reg8+reg17; double reg37=1.3660254037844385386*e.pos(3)[0]; reg14=reg15+reg14;
   double reg38=0.45534180126147951289*e.pos(3)[0]; double reg39=1.3660254037844385386*e.pos(2)[2]; double reg40=0.45534180126147951289*e.pos(0)[2]; double reg41=reg11+reg9; double reg42=0.12200846792814621817*e.pos(2)[1];
   double reg43=0.12200846792814621817*e.pos(2)[0]; double reg44=reg34+reg37; double reg45=0.45534180126147951289*e.pos(1)[2]; double reg46=reg27+reg26; double reg47=reg24+reg28;
   double reg48=1.3660254037844385386*e.pos(3)[2]; double reg49=0.12200846792814621817*e.pos(2)[2]; double reg50=1.3660254037844385386*e.pos(0)[2]; double reg51=0.36602540378443865451*e.pos(2)[1]; reg36=reg43+reg36;
   double reg52=0.36602540378443865451*e.pos(3)[0]; reg33=reg42+reg33; double reg53=0.36602540378443865451*e.pos(3)[1]; double reg54=reg40+reg31; double reg55=reg35+reg32;
   double reg56=reg30+reg29; double reg57=0.36602540378443865451*e.pos(2)[0]; double reg58=0.48803387171258487271*e.pos(4)[1]; reg12=reg12+reg25; reg14=reg14+reg38;
   double reg59=0.45534180126147951289*e.pos(3)[2]; reg41=reg41+reg39; double reg60=0.48803387171258487271*e.pos(4)[0]; double reg61=0.66666666666666670528*e.pos(5)[1]; double reg62=0.12200846792814621817*e.pos(3)[1];
   reg33=reg33+reg53; double reg63=reg58-reg46; double reg64=1.8213672050459180516*e.pos(4)[0]; reg36=reg36+reg52; reg12=reg12-reg58;
   double reg65=0.66666666666666670528*e.pos(5)[0]; double reg66=reg55+reg51; double reg67=reg45+reg50; double reg68=0.36602540378443865451*e.pos(2)[2]; double reg69=reg47+reg48;
   reg14=reg14-reg60; double reg70=reg60-reg44; double reg71=reg57+reg56; reg41=reg41+reg59; double reg72=0.48803387171258487271*e.pos(4)[2];
   double reg73=0.12200846792814621817*e.pos(3)[0]; double reg74=0.36602540378443865451*e.pos(3)[2]; reg54=reg49+reg54; double reg75=1.8213672050459180516*e.pos(4)[1]; double reg76=reg73+reg71;
   reg12=reg12+reg61; double reg77=reg72-reg69; reg70=reg65+reg70; double reg78=reg67+reg68; reg33=reg33-reg75;
   double reg79=1.8213672050459180516*e.pos(6)[1]; reg41=reg41-reg72; double reg80=0.66666666666666670528*e.pos(5)[2]; double reg81=reg62+reg66; double reg82=0.12200846792814621817*e.pos(3)[2];
   double reg83=1.8213672050459180516*e.pos(6)[0]; double reg84=1.8213672050459180516*e.pos(4)[2]; reg14=reg65+reg14; reg54=reg54+reg74; reg36=reg36-reg64;
   reg63=reg61+reg63; reg14=reg14-reg83; double reg85=reg5+reg8; reg77=reg80+reg77; double reg86=0.66666666666666670528*e.pos(7)[0];
   reg12=reg12-reg79; reg70=reg83+reg70; double reg87=0.66666666666666670528*e.pos(7)[1]; reg63=reg79+reg63; double reg88=reg7+reg10;
   double reg89=reg64-reg76; double reg90=0.48803387171258487271*e.pos(6)[1]; reg33=reg61+reg33; double reg91=1.8213672050459180516*e.pos(6)[2]; double reg92=reg82+reg78;
   double reg93=reg3+reg35; double reg94=reg1+reg30; double reg95=reg75-reg81; reg41=reg41+reg80; reg54=reg54-reg84;
   reg36=reg65+reg36; double reg96=0.48803387171258487271*e.pos(6)[0]; reg70=reg70-reg86; reg54=reg80+reg54; double reg97=0.48803387171258487271*e.pos(6)[2];
   reg17=reg4+reg17; reg20=reg6+reg20; reg89=reg65+reg89; reg12=reg12-reg87; reg33=reg33-reg90;
   reg36=reg36-reg96; reg63=reg63-reg87; reg14=reg14-reg86; reg95=reg61+reg95; reg4=reg84-reg92;
   reg6=reg18+reg40; reg77=reg91+reg77; reg42=reg42+reg88; reg43=reg43+reg85; reg15=reg15+reg94;
   reg41=reg41-reg91; reg13=reg13+reg93; reg61=reg9+reg45; reg65=0.66666666666666670528*e.pos(7)[2]; reg4=reg80+reg4;
   reg95=reg90+reg95; reg39=reg39+reg61; reg29=reg0+reg29; reg43=reg37+reg43; reg36=reg36-reg86;
   reg0=0.66666666666666670528*e.pos(4)[1]; reg62=reg13+reg62; reg77=reg77-reg65; reg49=reg49+reg6; reg33=reg33-reg87;
   reg13=0.66666666666666670528*e.pos(4)[0]; reg37=pow(reg14,2); reg42=reg26+reg42; reg41=reg41-reg65; reg26=pow(reg70,2);
   reg80=pow(reg12,2); reg31=reg16+reg31; reg32=reg2+reg32; reg20=reg21+reg20; reg17=reg23+reg17;
   reg89=reg96+reg89; reg15=reg73+reg15; reg54=reg54-reg97; reg2=pow(reg63,2); reg42=reg42-reg0;
   reg16=0.48803387171258487271*e.pos(5)[1]; reg73=pow(reg41,2); reg49=reg48+reg49; reg4=reg97+reg4; reg87=reg95-reg87;
   reg31=reg28+reg31; reg48=0.66666666666666670528*e.pos(4)[2]; reg82=reg39+reg82; reg15=reg15-reg13; reg39=reg53+reg20;
   reg95=1.8213672050459180516*e.pos(5)[0]; double reg98=1.8213672050459180516*e.pos(5)[1]; double reg99=pow(reg33,2); double reg100=pow(reg36,2); reg62=reg62-reg0;
   double reg101=reg52+reg17; reg54=reg54-reg65; reg86=reg89-reg86; reg80=reg37+reg80; reg32=reg51+reg32;
   reg2=reg26+reg2; reg50=reg11+reg50; reg11=pow(reg77,2); reg29=reg57+reg29; reg43=reg43-reg13;
   reg26=0.48803387171258487271*e.pos(5)[0]; reg37=1.8213672050459180516*e.pos(5)[2]; reg89=pow(reg87,2); reg82=reg82-reg48; reg43=reg43-reg26;
   reg65=reg4-reg65; reg4=0.66666666666666670528*e.pos(6)[1]; reg62=reg62-reg98; double reg102=reg38+reg29; reg11=reg2+reg11;
   reg2=0.66666666666666670528*e.pos(6)[0]; reg15=reg15-reg95; reg99=reg100+reg99; reg50=reg68+reg50; reg100=pow(reg54,2);
   double reg103=pow(reg86,2); reg42=reg42-reg16; double reg104=reg25+reg32; double reg105=0.48803387171258487271*e.pos(5)[2]; reg49=reg49-reg48;
   double reg106=reg13+reg101; reg73=reg80+reg73; reg80=reg0+reg39; double reg107=reg74+reg31; reg89=reg103+reg89;
   reg13=reg13+reg102; reg103=0.48803387171258487271*e.pos(7)[1]; reg73=pow(reg73,0.5); reg62=reg62+reg4; double reg108=reg48+reg107;
   double reg109=reg98-reg80; reg11=pow(reg11,0.5); double reg110=0.48803387171258487271*e.pos(7)[0]; reg15=reg15+reg2; double reg111=reg95-reg106;
   double reg112=reg59+reg50; reg100=reg99+reg100; reg0=reg0+reg104; reg82=reg82-reg37; reg99=0.66666666666666670528*e.pos(6)[2];
   double reg113=pow(reg65,2); reg42=reg4+reg42; double reg114=1.8213672050459180516*e.pos(7)[1]; reg43=reg43+reg2; double reg115=1.8213672050459180516*e.pos(7)[0];
   reg49=reg49-reg105; reg109=reg4+reg109; double reg116=reg63/reg11; double reg117=reg70/reg11; reg111=reg2+reg111;
   reg49=reg99+reg49; reg15=reg15-reg110; double reg118=reg12/reg73; double reg119=reg16-reg0; double reg120=1.8213672050459180516*e.pos(7)[2];
   reg100=pow(reg100,0.5); reg48=reg48+reg112; reg82=reg82+reg99; double reg121=0.48803387171258487271*e.pos(7)[2]; double reg122=reg26-reg13;
   reg42=reg42-reg114; reg43=reg43-reg115; double reg123=reg14/reg73; double reg124=reg37-reg108; reg62=reg62-reg103;
   reg113=reg89+reg113; reg124=reg99+reg124; reg122=reg2+reg122; reg2=reg116*reg42; reg111=reg110+reg111;
   reg119=reg4+reg119; reg4=reg117*reg43; reg89=reg33/reg100; double reg125=reg36/reg100; reg109=reg103+reg109;
   reg49=reg49-reg120; reg73=reg41/reg73; double reg126=reg118*reg62; reg113=pow(reg113,0.5); reg82=reg82-reg121;
   double reg127=reg123*reg15; double reg128=reg105-reg48; reg11=reg77/reg11; double reg129=reg73*reg82; reg2=reg4+reg2;
   reg4=reg125*reg111; reg126=reg127+reg126; reg124=reg121+reg124; reg127=reg86/reg113; double reg130=reg89*reg109;
   reg100=reg54/reg100; reg122=reg115+reg122; double reg131=reg87/reg113; reg119=reg114+reg119; double reg132=reg11*reg49;
   reg128=reg99+reg128; reg113=reg65/reg113; reg128=reg120+reg128; reg99=reg100*reg124; reg132=reg2+reg132;
   reg2=reg127*reg122; double reg133=reg131*reg119; reg129=reg126+reg129; reg130=reg4+reg130; reg99=reg130+reg99;
   reg4=reg117*reg132; reg126=reg116*reg132; reg130=reg118*reg129; double reg134=reg113*reg128; double reg135=reg123*reg129;
   reg133=reg2+reg133; reg2=reg125*reg99; double reg136=reg89*reg99; double reg137=reg73*reg129; reg130=reg62-reg130;
   reg134=reg133+reg134; reg133=reg11*reg132; reg135=reg15-reg135; reg126=reg42-reg126; reg4=reg43-reg4;
   double reg138=pow(reg135,2); reg137=reg82-reg137; double reg139=reg127*reg134; double reg140=reg100*reg99; reg136=reg109-reg136;
   double reg141=reg131*reg134; double reg142=pow(reg4,2); reg2=reg111-reg2; double reg143=pow(reg126,2); reg133=reg49-reg133;
   double reg144=pow(reg130,2); reg139=reg122-reg139; double reg145=reg113*reg134; reg141=reg119-reg141; double reg146=pow(reg2,2);
   double reg147=pow(reg136,2); reg140=reg124-reg140; double reg148=pow(reg137,2); reg144=reg138+reg144; reg138=pow(reg133,2);
   reg143=reg142+reg143; reg147=reg146+reg147; reg142=pow(reg140,2); reg138=reg143+reg138; reg148=reg144+reg148;
   reg145=reg128-reg145; reg143=pow(reg139,2); reg144=pow(reg141,2); reg148=pow(reg148,0.5); reg144=reg143+reg144;
   reg138=pow(reg138,0.5); reg142=reg147+reg142; reg143=pow(reg145,2); reg135=reg135/reg148; reg130=reg130/reg148;
   reg4=reg4/reg138; reg142=pow(reg142,0.5); reg126=reg126/reg138; reg143=reg144+reg143; reg43=reg43*reg4;
   reg42=reg42*reg126; reg138=reg133/reg138; reg116=reg63*reg116; reg117=reg70*reg117; reg133=reg12*reg130;
   reg144=reg14*reg135; reg148=reg137/reg148; reg130=reg62*reg130; reg135=reg15*reg135; reg118=reg12*reg118;
   reg123=reg14*reg123; reg143=pow(reg143,0.5); reg2=reg2/reg142; reg136=reg136/reg142; reg126=reg63*reg126;
   reg4=reg70*reg4; reg82=reg82*reg148; reg130=reg135+reg130; reg133=reg144+reg133; reg148=reg41*reg148;
   reg12=reg33*reg136; reg14=reg36*reg2; reg73=reg41*reg73; reg118=reg123+reg118; reg89=reg33*reg89;
   reg116=reg117+reg116; reg142=reg140/reg142; reg11=reg77*reg11; reg136=reg109*reg136; reg49=reg49*reg138;
   reg42=reg43+reg42; reg126=reg4+reg126; reg138=reg77*reg138; reg141=reg141/reg143; reg125=reg36*reg125;
   reg139=reg139/reg143; reg2=reg111*reg2; reg122=reg122*reg139; reg124=reg124*reg142; reg119=reg119*reg141;
   reg143=reg145/reg143; reg73=reg118+reg73; reg139=reg86*reg139; reg136=reg2+reg136; reg12=reg14+reg12;
   reg142=reg54*reg142; reg141=reg87*reg141; reg49=reg42+reg49; reg138=reg126+reg138; reg127=reg86*reg127;
   reg11=reg116+reg11; reg131=reg87*reg131; reg89=reg125+reg89; reg148=reg133+reg148; reg100=reg54*reg100;
   reg82=reg130+reg82; reg119=reg122+reg119; reg138=reg132*reg138; reg82=reg73*reg82; reg128=reg128*reg143;
   reg49=reg11*reg49; reg131=reg127+reg131; reg113=reg65*reg113; reg141=reg139+reg141; reg143=reg65*reg143;
   reg124=reg136+reg124; reg148=reg129*reg148; reg142=reg12+reg142; reg100=reg89+reg100; reg143=reg141+reg143;
   reg128=reg119+reg128; reg138=reg49-reg138; reg148=reg82-reg148; reg113=reg131+reg113; reg142=reg99*reg142;
   reg124=reg100*reg124; reg2=0.024056261216234395431*reg148; reg4=0.13144585576580215187*reg138; reg11=0.024056261216234395431*reg138; reg12=0.024056261216234409915*reg138;
   reg14=0.13144585576580215187*reg148; reg15=0.04166666666666666908*reg148; reg33=0.035220810900864524453*reg148; reg36=0.024056261216234409915*reg148; reg41=0.035220810900864524453*reg138;
   reg143=reg134*reg143; reg128=reg113*reg128; reg142=reg124-reg142; reg42=0.04166666666666666908*reg138; reg43=reg33+reg41;
   reg49=0.024056261216234409915*reg142; reg11=reg11-reg15; reg54=0.13144585576580215187*reg142; reg41=reg41+reg14; reg2=reg2-reg42;
   reg14=reg14+reg4; reg62=0.035220810900864524453*reg142; reg12=reg15+reg12; reg15=0.024056261216234395431*reg142; reg4=reg33+reg4;
   reg143=reg128-reg143; reg33=0.04166666666666666908*reg142; reg42=reg36+reg42; reg36=0.035220810900864524453*reg143; reg41=reg54+reg41;
   reg63=0.13144585576580215187*reg143; reg14=reg14+reg62; reg4=reg62+reg4; reg54=reg43+reg54; reg49=reg11-reg49;
   reg11=0.024056261216234409915*reg143; reg2=reg2-reg33; reg43=0.04166666666666666908*reg143; reg12=reg15-reg12; reg15=0.024056261216234395431*reg143;
   reg33=reg42+reg33; Ne(0,9)+=reg49-reg43; Ne(1,10)+=reg49-reg43; Ne(2,11)+=reg49-reg43; Ne(0,12)+=reg54+reg63;
   Ne(1,13)+=reg54+reg63; Ne(2,14)+=reg54+reg63; Ne(0,6)+=reg2-reg11; Ne(1,7)+=reg2-reg11; Ne(2,8)+=reg2-reg11;
   Ne(0,15)+=reg41+reg36; Ne(1,16)+=reg41+reg36; Ne(2,17)+=reg41+reg36; Ne(0,3)+=reg12-reg43; Ne(1,4)+=reg12-reg43;
   Ne(2,5)+=reg12-reg43; Ne(0,0)+=reg15-reg33; Ne(1,1)+=reg15-reg33; Ne(2,2)+=reg15-reg33; Ne(0,18)+=reg36+reg14;
   Ne(1,19)+=reg36+reg14; Ne(2,20)+=reg36+reg14; Ne(0,21)+=reg63+reg4; Ne(1,22)+=reg63+reg4; Ne(2,23)+=reg63+reg4;

}

};
} // namespace LMT

