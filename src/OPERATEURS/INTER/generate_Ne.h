
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
   
   double reg0=0.56758792732771999991*e.pos(1)[0]; double reg1=0.78379396366385999995*e.pos(1)[0]; double reg2=0.56758792732771999991*e.pos(0)[0]; double reg3=0.56758792732771999991*e.pos(1)[1]; double reg4=0.56758792732771999991*e.pos(0)[1];
   double reg5=0.78379396366385999995*e.pos(1)[1]; double reg6=0.78379396366385999995*e.pos(0)[0]; double reg7=0.78379396366385999995*e.pos(0)[1]; double reg8=1.3513818909915799999*e.pos(3)[0]; double reg9=reg5+reg4;
   double reg10=reg1+reg2; double reg11=1.3513818909915799999*e.pos(3)[1]; double reg12=reg3+reg7; double reg13=reg0+reg6; double reg14=0.56758792732771999991*e.pos(1)[2];
   double reg15=0.78379396366385999995*e.pos(0)[2]; double reg16=0.78379396366385999995*e.pos(1)[2]; double reg17=0.56758792732771999991*e.pos(0)[2]; reg9=reg9-reg11; double reg18=1.78379396366386*e.pos(4)[1];
   double reg19=reg8-reg13; double reg20=2.2673902919218320001*e.pos(0)[0]; double reg21=reg16+reg17; double reg22=1.78379396366386*e.pos(4)[0]; reg10=reg10-reg8;
   double reg23=reg11-reg12; double reg24=0.63369514596091600003*e.pos(1)[0]; double reg25=0.63369514596091600003*e.pos(1)[1]; double reg26=reg14+reg15; double reg27=1.3513818909915799999*e.pos(3)[2];
   double reg28=2.2673902919218320001*e.pos(0)[1]; double reg29=0.63369514596091600003*e.pos(0)[1]; double reg30=2.2673902919218320001*e.pos(1)[0]; reg5=reg5-reg7; double reg31=2.2673902919218320001*e.pos(1)[1];
   double reg32=reg24+reg20; double reg33=0.63369514596091600003*e.pos(0)[0]; double reg34=2.9010854378827480001*e.pos(3)[0]; double reg35=2.9010854378827480001*e.pos(3)[1]; double reg36=2.2673902919218320001*e.pos(0)[2];
   double reg37=0.63369514596091600003*e.pos(1)[2]; double reg38=reg25+reg28; double reg39=0.43241207267228000009*e.pos(4)[0]; reg23=reg18+reg23; reg1=reg1-reg6;
   double reg40=reg27-reg26; double reg41=0.43241207267228000009*e.pos(4)[1]; double reg42=1.78379396366386*e.pos(4)[2]; reg21=reg21-reg27; reg19=reg22+reg19;
   reg10=reg10+reg22; double reg43=1.78379396366386*e.pos(5)[1]; reg9=reg9+reg18; double reg44=1.78379396366386*e.pos(5)[0]; double reg45=0.366304854039084*e.pos(4)[0];
   reg23=reg23-reg43; reg40=reg42+reg40; double reg46=reg35-reg38; double reg47=reg34-reg32; reg19=reg19-reg44;
   reg1=reg39+reg1; double reg48=0.43241207267228000009*e.pos(4)[2]; reg16=reg16-reg15; double reg49=0.43241207267228000009*e.pos(5)[1]; reg5=reg41+reg5;
   double reg50=0.43241207267228000009*e.pos(5)[0]; reg44=reg10-reg44; reg43=reg9-reg43; reg21=reg21+reg42; reg9=1.78379396366386*e.pos(5)[2];
   reg10=0.63369514596091600003*e.pos(0)[2]; double reg51=2.2673902919218320001*e.pos(1)[2]; reg31=reg31+reg29; reg30=reg30+reg33; double reg52=2.9010854378827480001*e.pos(3)[2];
   double reg53=reg37+reg36; double reg54=0.366304854039084*e.pos(4)[1]; reg51=reg51+reg10; reg50=reg1-reg50; reg1=reg33-reg24;
   double reg55=0.366304854039084*e.pos(5)[0]; reg21=reg21-reg9; double reg56=2.710505431213761085e-20*e.pos(3)[0]; double reg57=pow(reg43,2); double reg58=pow(reg19,2);
   reg9=reg40-reg9; reg40=pow(reg44,2); double reg59=2.710505431213761085e-20*e.pos(3)[1]; double reg60=reg29-reg25; double reg61=pow(reg23,2);
   double reg62=0.366304854039084*e.pos(5)[1]; reg46=reg46+reg54; double reg63=0.78379396366385999995*e.pos(2)[1]; reg30=reg30-reg34; double reg64=0.78379396366385999995*e.pos(2)[0];
   double reg65=0.366304854039084*e.pos(4)[2]; reg31=reg31-reg35; double reg66=0.43241207267228000009*e.pos(5)[2]; reg16=reg48+reg16; reg47=reg47+reg45;
   reg49=reg5-reg49; reg5=reg52-reg53; double reg67=0.366304854039084*e.pos(5)[2]; double reg68=2.710505431213761085e-20*e.pos(3)[2]; double reg69=pow(reg50,2);
   reg5=reg5+reg65; double reg70=reg10-reg37; reg60=reg60-reg59; reg30=reg45+reg30; reg51=reg51-reg52;
   double reg71=3.2673902919218320001*e.pos(4)[0]; reg1=reg1-reg56; reg46=reg46-reg62; double reg72=3.2673902919218320001*e.pos(4)[1]; double reg73=0.56758792732771999991*e.pos(2)[1];
   double reg74=0.56758792732771999991*e.pos(2)[0]; reg66=reg16-reg66; reg16=pow(reg49,2); reg47=reg47-reg55; reg31=reg54+reg31;
   reg61=reg58+reg61; reg58=reg2+reg64; double reg75=1.78379396366386*e.pos(3)[0]; double reg76=pow(reg9,2); double reg77=pow(reg21,2);
   reg57=reg40+reg57; reg40=reg4+reg63; double reg78=1.78379396366386*e.pos(3)[1]; double reg79=0.78379396366385999995*e.pos(2)[2]; reg1=reg71+reg1;
   double reg80=3.2673902919218320001*e.pos(5)[0]; reg60=reg72+reg60; reg77=reg57+reg77; reg57=3.2673902919218320001*e.pos(4)[2]; reg16=reg69+reg16;
   reg69=pow(reg66,2); double reg81=3.2673902919218320001*e.pos(5)[1]; reg51=reg65+reg51; reg70=reg70-reg68; reg62=reg31-reg62;
   reg31=reg6+reg74; double reg82=reg7+reg73; reg55=reg30-reg55; reg58=reg58-reg75; reg30=0.56758792732771999991*e.pos(2)[2];
   double reg83=0.63369514596091600003*e.pos(2)[1]; double reg84=0.63369514596091600003*e.pos(2)[0]; reg5=reg5-reg67; reg40=reg40-reg78; double reg85=pow(reg46,2);
   double reg86=reg17+reg79; double reg87=1.78379396366386*e.pos(3)[2]; double reg88=pow(reg47,2); reg76=reg61+reg76; reg61=0.43241207267228000009*e.pos(3)[1];
   reg7=reg63-reg7; reg6=reg64-reg6; reg63=0.43241207267228000009*e.pos(3)[0]; reg79=reg79-reg15; reg67=reg51-reg67;
   reg69=reg16+reg69; reg16=pow(reg62,2); reg61=reg7-reg61; reg75=reg75+reg31; reg7=pow(reg55,2);
   reg78=reg78+reg82; reg63=reg6-reg63; reg15=reg15+reg30; reg58=reg22+reg58; reg6=1.3513818909915799999*e.pos(5)[0];
   reg51=0.63369514596091600003*e.pos(2)[2]; reg64=0.366304854039084*e.pos(3)[1]; reg28=reg28+reg83; double reg89=0.366304854039084*e.pos(3)[0]; reg20=reg20+reg84;
   double reg90=pow(reg5,2); reg40=reg18+reg40; double reg91=1.3513818909915799999*e.pos(5)[1]; reg85=reg88+reg85; reg76=pow(reg76,0.5);
   reg86=reg86-reg87; reg77=pow(reg77,0.5); reg88=0.43241207267228000009*e.pos(3)[2]; reg80=reg1-reg80; reg81=reg60-reg81;
   reg70=reg57+reg70; reg1=3.2673902919218320001*e.pos(5)[2]; reg18=reg18-reg78; reg39=reg63+reg39; reg87=reg87+reg15;
   reg60=1.3513818909915799999*e.pos(5)[2]; reg88=reg79-reg88; reg58=reg58-reg6; reg63=0.366304854039084*e.pos(3)[2]; reg36=reg36+reg51;
   reg79=reg28+reg64; double reg92=2.2673902919218320001*e.pos(2)[0]; double reg93=reg20+reg89; reg86=reg42+reg86; double reg94=pow(reg80,2);
   double reg95=2.2673902919218320001*e.pos(2)[1]; reg90=reg85+reg90; reg1=reg70-reg1; reg40=reg40-reg91; reg70=pow(reg81,2);
   reg85=reg23/reg76; double reg96=reg19/reg76; reg41=reg61+reg41; reg69=pow(reg69,0.5); reg61=reg29-reg83;
   double reg97=3.2673902919218320001*e.pos(3)[1]; reg22=reg22-reg75; double reg98=pow(reg67,2); reg16=reg7+reg16; reg7=reg44/reg77;
   double reg99=reg33-reg84; double reg100=3.2673902919218320001*e.pos(3)[0]; double reg101=reg43/reg77; double reg102=reg49/reg69; double reg103=reg85*reg41;
   double reg104=reg50/reg69; reg90=pow(reg90,0.5); double reg105=pow(reg1,2); double reg106=reg45-reg93; reg95=reg29+reg95;
   double reg107=reg7*reg58; reg76=reg9/reg76; double reg108=2.2673902919218320001*e.pos(2)[2]; double reg109=reg101*reg40; reg98=reg16+reg98;
   reg70=reg94+reg70; reg86=reg86-reg60; reg100=reg99-reg100; reg16=reg96*reg39; reg18=reg91+reg18;
   reg22=reg6+reg22; reg94=reg10-reg51; reg99=3.2673902919218320001*e.pos(3)[2]; reg42=reg42-reg87; reg97=reg61-reg97;
   reg61=reg36+reg63; reg77=reg21/reg77; reg92=reg33+reg92; double reg110=2.9010854378827480001*e.pos(5)[0]; double reg111=2.9010854378827480001*e.pos(5)[1];
   reg48=reg88+reg48; reg88=reg54-reg79; double reg112=reg77*reg86; double reg113=5.42101086242752217e-20*e.pos(5)[1]; reg72=reg97+reg72;
   reg106=reg106+reg110; reg105=reg70+reg105; reg98=pow(reg98,0.5); reg108=reg10+reg108; reg109=reg107+reg109;
   reg70=5.42101086242752217e-20*e.pos(5)[0]; reg99=reg94-reg99; reg71=reg100+reg71; reg89=reg92-reg89; reg92=reg46/reg90;
   reg94=2.9010854378827480001*e.pos(5)[2]; reg64=reg95-reg64; reg95=reg65-reg61; reg97=reg47/reg90; reg88=reg88+reg111;
   reg103=reg16+reg103; reg16=reg76*reg48; reg69=reg66/reg69; reg100=reg104*reg22; reg107=reg102*reg18;
   reg42=reg60+reg42; reg63=reg108-reg63; reg90=reg5/reg90; reg64=reg54+reg64; reg89=reg45+reg89;
   reg45=reg97*reg106; reg54=reg92*reg88; reg95=reg95+reg94; reg105=pow(reg105,0.5); reg16=reg103+reg16;
   reg103=reg69*reg42; reg107=reg100+reg107; reg100=reg55/reg98; reg108=reg62/reg98; double reg114=5.42101086242752217e-20*e.pos(5)[2];
   reg57=reg99+reg57; reg72=reg72-reg113; reg71=reg71-reg70; reg112=reg109+reg112; reg99=reg96*reg16;
   reg103=reg107+reg103; reg107=reg80/reg105; reg109=reg81/reg105; double reg115=reg85*reg16; reg54=reg45+reg54;
   reg45=reg101*reg112; reg89=reg89-reg110; reg98=reg67/reg98; reg64=reg64-reg111; reg57=reg57-reg114;
   double reg116=reg90*reg95; double reg117=reg108*reg72; double reg118=reg7*reg112; double reg119=reg100*reg71; reg63=reg65+reg63;
   reg116=reg54+reg116; reg105=reg1/reg105; reg54=reg104*reg103; reg99=reg39-reg99; reg65=reg77*reg112;
   double reg120=reg76*reg16; reg115=reg41-reg115; reg45=reg40-reg45; double reg121=reg98*reg57; reg118=reg58-reg118;
   reg117=reg119+reg117; reg63=reg63-reg94; reg119=reg102*reg103; double reg122=reg107*reg89; double reg123=reg109*reg64;
   double reg124=pow(reg99,2); reg120=reg48-reg120; double reg125=reg69*reg103; double reg126=pow(reg115,2); reg65=reg86-reg65;
   double reg127=pow(reg45,2); reg121=reg117+reg121; reg117=reg92*reg116; double reg128=pow(reg118,2); reg123=reg122+reg123;
   reg122=reg97*reg116; reg119=reg18-reg119; double reg129=reg105*reg63; reg54=reg22-reg54; double reg130=reg100*reg121;
   double reg131=reg108*reg121; reg127=reg128+reg127; reg128=pow(reg65,2); reg124=reg126+reg124; reg129=reg123+reg129;
   reg123=pow(reg120,2); reg126=pow(reg54,2); reg125=reg42-reg125; double reg132=pow(reg119,2); reg122=reg106-reg122;
   reg117=reg88-reg117; double reg133=reg90*reg116; double reg134=pow(reg122,2); double reg135=reg107*reg129; double reg136=pow(reg117,2);
   reg130=reg71-reg130; reg133=reg95-reg133; reg131=reg72-reg131; double reg137=reg98*reg121; double reg138=reg109*reg129;
   double reg139=pow(reg125,2); reg123=reg124+reg123; reg132=reg126+reg132; reg128=reg127+reg128; reg124=pow(reg130,2);
   reg126=pow(reg131,2); reg139=reg132+reg139; reg137=reg57-reg137; reg135=reg89-reg135; reg136=reg134+reg136;
   reg127=pow(reg133,2); reg123=pow(reg123,0.5); reg138=reg64-reg138; reg132=reg105*reg129; reg128=pow(reg128,0.5);
   reg99=reg99/reg123; reg115=reg115/reg123; reg134=pow(reg135,2); reg118=reg118/reg128; double reg140=pow(reg138,2);
   reg126=reg124+reg126; reg124=pow(reg137,2); reg132=reg63-reg132; reg139=pow(reg139,0.5); reg45=reg45/reg128;
   reg127=reg136+reg127; reg140=reg134+reg140; reg96=reg19*reg96; reg85=reg23*reg85; reg123=reg120/reg123;
   reg19=reg19*reg99; reg23=reg23*reg115; reg99=reg39*reg99; reg124=reg126+reg124; reg39=pow(reg132,2);
   reg115=reg41*reg115; reg40=reg40*reg45; reg127=pow(reg127,0.5); reg101=reg43*reg101; reg7=reg44*reg7;
   reg128=reg65/reg128; reg54=reg54/reg139; reg119=reg119/reg139; reg44=reg44*reg118; reg45=reg43*reg45;
   reg118=reg58*reg118; reg117=reg117/reg127; reg122=reg122/reg127; reg115=reg99+reg115; reg124=pow(reg124,0.5);
   reg101=reg7+reg101; reg77=reg21*reg77; reg102=reg49*reg102; reg49=reg49*reg119; reg104=reg50*reg104;
   reg40=reg118+reg40; reg50=reg50*reg54; reg86=reg86*reg128; reg48=reg48*reg123; reg54=reg22*reg54;
   reg119=reg18*reg119; reg139=reg125/reg139; reg39=reg140+reg39; reg23=reg19+reg23; reg123=reg9*reg123;
   reg76=reg9*reg76; reg85=reg96+reg85; reg128=reg21*reg128; reg45=reg44+reg45; reg48=reg115+reg48;
   reg69=reg66*reg69; reg102=reg104+reg102; reg77=reg101+reg77; reg7=reg46*reg117; reg9=reg47*reg122;
   reg127=reg133/reg127; reg117=reg88*reg117; reg122=reg106*reg122; reg66=reg66*reg139; reg49=reg50+reg49;
   reg97=reg47*reg97; reg92=reg46*reg92; reg119=reg54+reg119; reg76=reg85+reg76; reg139=reg42*reg139;
   reg128=reg45+reg128; reg123=reg23+reg123; reg86=reg40+reg86; reg130=reg130/reg124; reg131=reg131/reg124;
   reg39=pow(reg39,0.5); reg95=reg95*reg127; reg124=reg137/reg124; reg117=reg122+reg117; reg72=reg72*reg131;
   reg48=reg76*reg48; reg71=reg71*reg130; reg138=reg138/reg39; reg100=reg55*reg100; reg108=reg62*reg108;
   reg86=reg77*reg86; reg135=reg135/reg39; reg123=reg16*reg123; reg139=reg119+reg139; reg90=reg5*reg90;
   reg92=reg97+reg92; reg66=reg49+reg66; reg128=reg112*reg128; reg69=reg102+reg69; reg130=reg55*reg130;
   reg131=reg62*reg131; reg7=reg9+reg7; reg127=reg5*reg127; reg5=reg80*reg135; reg123=reg48-reg123;
   reg90=reg92+reg90; reg131=reg130+reg131; reg139=reg69*reg139; reg128=reg86-reg128; reg9=reg67*reg124;
   reg64=reg64*reg138; reg66=reg103*reg66; reg95=reg117+reg95; reg124=reg57*reg124; reg39=reg132/reg39;
   reg72=reg71+reg72; reg107=reg80*reg107; reg109=reg81*reg109; reg127=reg7+reg127; reg135=reg89*reg135;
   reg108=reg100+reg108; reg98=reg67*reg98; reg138=reg81*reg138; reg7=reg1*reg39; reg16=0.005384432036113586778*reg128;
   reg18=0.088847818743090689935*reg128; reg19=0.088847818743090689935*reg123; reg138=reg5+reg138; reg64=reg135+reg64; reg39=reg63*reg39;
   reg5=0.021537728144454347112*reg123; reg21=0.021537728144454347112*reg128; reg22=0.009463616120767210603*reg123; reg105=reg1*reg105; reg109=reg107+reg109;
   reg1=0.009463616120767210603*reg128; reg127=reg116*reg127; reg95=reg90*reg95; reg66=reg139-reg66; reg9=reg131+reg9;
   reg23=0.005384432036113586778*reg123; reg124=reg72+reg124; reg98=reg108+reg98; reg19=reg21+reg19; reg22=reg16+reg22;
   reg16=reg23+reg16; reg40=0.009463616120767210603*reg66; reg41=0.021537728144454347112*reg66; reg18=reg5+reg18; reg42=0.005384432036113586778*reg66;
   reg5=reg21+reg5; reg23=reg1+reg23; reg1=0.088847818743090689935*reg66; reg105=reg109+reg105; reg7=reg138+reg7;
   reg127=reg95-reg127; reg9=reg121*reg9; reg124=reg98*reg124; reg39=reg64+reg39; reg21=0.016449618187943419918*reg127;
   reg1=reg5+reg1; reg9=reg124-reg9; reg23=reg23+reg42; reg39=reg105*reg39; reg40=reg16+reg40;
   reg18=reg18+reg41; reg5=0.028457289286966203713*reg127; reg16=0.0018441552587796664112*reg127; reg43=0.0041124045469858549794*reg127; reg22=reg42+reg22;
   reg7=reg129*reg7; reg19=reg41+reg19; reg19=reg21+reg19; reg16=reg18+reg16; reg18=0.0018441552587796664109*reg9;
   reg41=0.016449618187943419918*reg9; reg42=0.016449618187943419916*reg9; reg44=0.004112404546985854979*reg9; reg23=reg5-reg23; reg7=reg39-reg7;
   reg5=0.028457289286966203713*reg9; reg22=reg22+reg43; reg40=reg43+reg40; reg39=0.0041124045469858549794*reg9; reg21=reg1+reg21;
   reg1=0.028457289286966203713*reg7; reg39=reg40+reg39; reg22=reg5-reg22; reg18=reg19+reg18; reg5=0.0018441552587796664111*reg7;
   reg41=reg16+reg41; reg16=0.016449618187943419918*reg7; reg44=reg23-reg44; reg19=0.0041124045469858549794*reg7; reg42=reg21+reg42;
   Ne(0,15)+=reg16+reg18; Ne(1,16)+=reg16+reg18; Ne(2,17)+=reg16+reg18; Ne(0,6)+=reg1-reg39; Ne(1,7)+=reg1-reg39;
   Ne(2,8)+=reg1-reg39; Ne(0,0)+=reg44-reg19; Ne(1,1)+=reg44-reg19; Ne(2,2)+=reg44-reg19; Ne(0,12)+=reg41+reg16;
   Ne(1,13)+=reg41+reg16; Ne(2,14)+=reg41+reg16; Ne(0,9)+=reg42+reg5; Ne(1,10)+=reg42+reg5; Ne(2,11)+=reg42+reg5;
   Ne(0,3)+=reg22-reg19; Ne(1,4)+=reg22-reg19; Ne(2,5)+=reg22-reg19;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.12200846792814621817*e.pos(1)[1]; double reg1=0.36602540378443865451*e.pos(0)[1]; double reg2=0.36602540378443865451*e.pos(0)[0]; double reg3=0.12200846792814621817*e.pos(1)[0]; double reg4=0.12200846792814621817*e.pos(0)[1];
   double reg5=0.36602540378443865451*e.pos(1)[1]; double reg6=0.12200846792814621817*e.pos(0)[0]; double reg7=0.36602540378443865451*e.pos(1)[0]; double reg8=1.3660254037844385386*e.pos(1)[0]; double reg9=0.36602540378443865451*e.pos(1)[2];
   double reg10=0.12200846792814621817*e.pos(0)[2]; double reg11=reg3+reg2; double reg12=0.45534180126147951289*e.pos(2)[0]; double reg13=1.3660254037844385386*e.pos(2)[0]; double reg14=0.45534180126147951289*e.pos(0)[0];
   double reg15=0.45534180126147951289*e.pos(2)[1]; double reg16=reg5+reg4; double reg17=1.3660254037844385386*e.pos(1)[1]; double reg18=reg7+reg6; double reg19=0.12200846792814621817*e.pos(1)[2];
   double reg20=0.36602540378443865451*e.pos(0)[2]; double reg21=1.3660254037844385386*e.pos(2)[1]; double reg22=reg0+reg1; double reg23=0.45534180126147951289*e.pos(0)[1]; double reg24=0.45534180126147951289*e.pos(3)[1];
   reg22=reg22+reg21; double reg25=0.45534180126147951289*e.pos(1)[1]; double reg26=0.45534180126147951289*e.pos(1)[0]; double reg27=0.45534180126147951289*e.pos(2)[2]; double reg28=reg9+reg10;
   double reg29=1.3660254037844385386*e.pos(3)[1]; double reg30=reg16+reg15; double reg31=reg19+reg20; double reg32=1.3660254037844385386*e.pos(2)[2]; double reg33=1.3660254037844385386*e.pos(3)[0];
   double reg34=reg12+reg18; double reg35=1.3660254037844385386*e.pos(0)[0]; double reg36=0.12200846792814621817*e.pos(2)[0]; double reg37=0.12200846792814621817*e.pos(2)[1]; double reg38=0.45534180126147951289*e.pos(0)[2];
   double reg39=0.45534180126147951289*e.pos(3)[0]; reg11=reg13+reg11; double reg40=1.3660254037844385386*e.pos(0)[1]; double reg41=reg14+reg8; double reg42=1.3660254037844385386*e.pos(1)[2];
   double reg43=reg23+reg17; double reg44=reg34+reg33; double reg45=reg38+reg42; double reg46=reg25+reg40; double reg47=reg26+reg35;
   double reg48=0.45534180126147951289*e.pos(1)[2]; double reg49=0.36602540378443865451*e.pos(2)[0]; double reg50=0.36602540378443865451*e.pos(3)[1]; reg43=reg37+reg43; double reg51=0.36602540378443865451*e.pos(3)[0];
   reg41=reg36+reg41; double reg52=0.36602540378443865451*e.pos(2)[1]; double reg53=reg30+reg29; double reg54=1.3660254037844385386*e.pos(0)[2]; double reg55=0.12200846792814621817*e.pos(2)[2];
   double reg56=1.3660254037844385386*e.pos(3)[2]; double reg57=reg28+reg27; reg11=reg11+reg39; double reg58=0.48803387171258487271*e.pos(4)[0]; double reg59=0.45534180126147951289*e.pos(3)[2];
   reg22=reg22+reg24; double reg60=0.48803387171258487271*e.pos(4)[1]; reg31=reg31+reg32; double reg61=1.8213672050459180516*e.pos(4)[1]; reg43=reg43+reg50;
   double reg62=1.8213672050459180516*e.pos(4)[0]; reg41=reg41+reg51; reg31=reg31+reg59; double reg63=0.66666666666666670528*e.pos(5)[0]; reg22=reg22-reg60;
   double reg64=reg46+reg52; double reg65=0.66666666666666670528*e.pos(5)[1]; double reg66=reg60-reg53; double reg67=reg57+reg56; double reg68=reg48+reg54;
   double reg69=0.36602540378443865451*e.pos(2)[2]; reg11=reg11-reg58; double reg70=reg58-reg44; double reg71=0.48803387171258487271*e.pos(4)[2]; double reg72=0.12200846792814621817*e.pos(3)[0];
   double reg73=0.36602540378443865451*e.pos(3)[2]; reg45=reg55+reg45; double reg74=0.12200846792814621817*e.pos(3)[1]; double reg75=reg49+reg47; reg45=reg45+reg73;
   double reg76=1.8213672050459180516*e.pos(4)[2]; double reg77=reg74+reg64; double reg78=reg68+reg69; reg11=reg63+reg11; double reg79=1.8213672050459180516*e.pos(6)[0];
   double reg80=0.66666666666666670528*e.pos(5)[2]; reg31=reg31-reg71; reg70=reg63+reg70; double reg81=reg71-reg67; double reg82=1.8213672050459180516*e.pos(6)[1];
   reg22=reg22+reg65; reg66=reg65+reg66; double reg83=reg72+reg75; double reg84=0.12200846792814621817*e.pos(3)[2]; reg43=reg43-reg61;
   reg41=reg41-reg62; reg43=reg65+reg43; double reg85=reg62-reg83; reg22=reg22-reg82; double reg86=0.66666666666666670528*e.pos(7)[1];
   double reg87=reg1+reg25; reg81=reg80+reg81; reg70=reg79+reg70; double reg88=1.8213672050459180516*e.pos(6)[2]; reg31=reg31+reg80;
   double reg89=0.66666666666666670528*e.pos(7)[0]; reg11=reg11-reg79; double reg90=0.48803387171258487271*e.pos(6)[1]; double reg91=reg7+reg14; reg66=reg82+reg66;
   double reg92=reg5+reg23; double reg93=reg61-reg77; reg41=reg63+reg41; double reg94=0.48803387171258487271*e.pos(6)[0]; reg45=reg45-reg76;
   double reg95=reg84+reg78; double reg96=reg2+reg26; reg17=reg4+reg17; reg31=reg31-reg88; reg93=reg65+reg93;
   reg70=reg70-reg89; reg41=reg41-reg94; reg85=reg63+reg85; reg22=reg22-reg86; reg4=0.48803387171258487271*e.pos(6)[2];
   reg45=reg80+reg45; reg66=reg66-reg86; reg11=reg11-reg89; reg43=reg43-reg90; reg21=reg21+reg87;
   reg63=reg9+reg38; reg65=0.66666666666666670528*e.pos(7)[2]; reg37=reg37+reg92; double reg97=reg76-reg95; reg36=reg36+reg91;
   reg13=reg13+reg96; reg81=reg88+reg81; double reg98=reg20+reg48; reg8=reg6+reg8; reg8=reg12+reg8;
   reg17=reg15+reg17; reg43=reg43-reg86; reg55=reg55+reg63; reg6=pow(reg11,2); reg37=reg29+reg37;
   reg42=reg10+reg42; reg97=reg80+reg97; reg85=reg94+reg85; reg36=reg33+reg36; reg93=reg90+reg93;
   reg45=reg45-reg4; reg35=reg3+reg35; reg81=reg81-reg65; reg41=reg41-reg89; reg40=reg0+reg40;
   reg32=reg32+reg98; reg0=0.66666666666666670528*e.pos(4)[1]; reg74=reg21+reg74; reg3=0.66666666666666670528*e.pos(4)[0]; reg72=reg13+reg72;
   reg10=pow(reg70,2); reg13=pow(reg22,2); reg31=reg31-reg65; reg21=pow(reg66,2); reg37=reg37-reg0;
   reg29=0.48803387171258487271*e.pos(5)[1]; reg54=reg19+reg54; reg55=reg56+reg55; reg97=reg4+reg97; reg86=reg93-reg86;
   reg19=0.66666666666666670528*e.pos(4)[2]; reg33=pow(reg41,2); reg84=reg32+reg84; reg32=1.8213672050459180516*e.pos(5)[1]; reg74=reg74-reg0;
   reg56=pow(reg43,2); reg80=pow(reg31,2); reg93=1.8213672050459180516*e.pos(5)[0]; reg72=reg72-reg3; reg45=reg45-reg65;
   reg42=reg27+reg42; reg89=reg85-reg89; reg85=reg51+reg8; double reg99=reg50+reg17; reg35=reg49+reg35;
   double reg100=pow(reg81,2); reg40=reg52+reg40; reg36=reg36-reg3; double reg101=0.48803387171258487271*e.pos(5)[0]; reg13=reg6+reg13;
   reg21=reg10+reg21; reg36=reg36-reg101; reg65=reg97-reg65; reg6=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg19;
   reg10=reg24+reg40; reg97=0.66666666666666670528*e.pos(6)[1]; reg74=reg74-reg32; double reg102=pow(reg86,2); reg56=reg33+reg56;
   reg33=0.66666666666666670528*e.pos(6)[0]; reg72=reg72-reg93; double reg103=reg39+reg35; double reg104=pow(reg45,2); double reg105=0.48803387171258487271*e.pos(5)[2];
   reg55=reg55-reg19; double reg106=pow(reg89,2); reg54=reg69+reg54; reg80=reg13+reg80; reg100=reg21+reg100;
   reg13=reg3+reg85; reg21=reg0+reg99; double reg107=reg73+reg42; reg37=reg37-reg29; reg102=reg106+reg102;
   reg106=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg6; reg80=pow(reg80,0.5); double reg108=reg19+reg107; reg3=reg3+reg103;
   double reg109=reg32-reg21; double reg110=0.48803387171258487271*e.pos(7)[1]; reg74=reg74+reg97; double reg111=reg93-reg13; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg72=reg72+reg33; reg100=pow(reg100,0.5); reg104=reg56+reg104; reg0=reg0+reg10; reg37=reg97+reg37;
   reg56=1.8213672050459180516*e.pos(7)[1]; double reg113=1.8213672050459180516*e.pos(7)[0]; double reg114=reg59+reg54; reg36=reg33+reg36; reg55=reg55-reg105;
   double reg115=pow(reg65,2); double reg116=reg6-reg108; double reg117=reg70/reg100; reg104=pow(reg104,0.5); reg55=reg106+reg55;
   double reg118=1.8213672050459180516*e.pos(7)[2]; double reg119=reg66/reg100; reg111=reg33+reg111; reg72=reg72-reg112; reg36=reg36-reg113;
   double reg120=reg101-reg3; double reg121=reg11/reg80; reg74=reg74-reg110; reg37=reg37-reg56; reg109=reg97+reg109;
   double reg122=reg22/reg80; reg115=reg102+reg115; reg102=reg29-reg0; reg19=reg19+reg114; reg84=reg84+reg106;
   double reg123=0.48803387171258487271*e.pos(7)[2]; reg116=reg106+reg116; reg80=reg31/reg80; double reg124=reg43/reg104; double reg125=reg41/reg104;
   reg102=reg97+reg102; reg55=reg55-reg118; reg97=reg119*reg37; reg84=reg84-reg123; reg111=reg112+reg111;
   double reg126=reg121*reg72; reg120=reg33+reg120; reg100=reg81/reg100; reg33=reg117*reg36; reg115=pow(reg115,0.5);
   double reg127=reg105-reg19; double reg128=reg122*reg74; reg109=reg110+reg109; reg102=reg56+reg102; reg128=reg126+reg128;
   reg104=reg45/reg104; reg127=reg106+reg127; reg116=reg123+reg116; reg106=reg89/reg115; reg97=reg33+reg97;
   reg33=reg124*reg109; reg126=reg80*reg84; reg120=reg113+reg120; double reg129=reg100*reg55; double reg130=reg86/reg115;
   double reg131=reg125*reg111; reg33=reg131+reg33; reg127=reg118+reg127; reg131=reg104*reg116; reg129=reg97+reg129;
   reg126=reg128+reg126; reg97=reg130*reg102; reg115=reg65/reg115; reg128=reg106*reg120; double reg132=reg121*reg126;
   double reg133=reg122*reg126; double reg134=reg115*reg127; double reg135=reg119*reg129; reg131=reg33+reg131; reg33=reg117*reg129;
   reg97=reg128+reg97; reg128=reg124*reg131; double reg136=reg125*reg131; double reg137=reg80*reg126; reg132=reg72-reg132;
   reg134=reg97+reg134; reg97=reg100*reg129; reg133=reg74-reg133; reg33=reg36-reg33; reg135=reg37-reg135;
   double reg138=pow(reg132,2); reg97=reg55-reg97; double reg139=reg104*reg131; double reg140=pow(reg135,2); reg128=reg109-reg128;
   reg136=reg111-reg136; double reg141=pow(reg33,2); double reg142=pow(reg133,2); double reg143=reg106*reg134; double reg144=reg130*reg134;
   reg137=reg84-reg137; reg144=reg102-reg144; reg143=reg120-reg143; double reg145=reg115*reg134; double reg146=pow(reg136,2);
   double reg147=pow(reg128,2); reg139=reg116-reg139; double reg148=pow(reg137,2); reg142=reg138+reg142; reg138=pow(reg97,2);
   reg140=reg141+reg140; reg147=reg146+reg147; reg141=pow(reg139,2); reg138=reg140+reg138; reg148=reg142+reg148;
   reg145=reg127-reg145; reg140=pow(reg143,2); reg142=pow(reg144,2); reg148=pow(reg148,0.5); reg142=reg140+reg142;
   reg138=pow(reg138,0.5); reg141=reg147+reg141; reg140=pow(reg145,2); reg132=reg132/reg148; reg133=reg133/reg148;
   reg33=reg33/reg138; reg141=pow(reg141,0.5); reg135=reg135/reg138; reg140=reg142+reg140; reg36=reg36*reg33;
   reg37=reg37*reg135; reg138=reg97/reg138; reg119=reg66*reg119; reg117=reg70*reg117; reg97=reg22*reg133;
   reg142=reg11*reg132; reg148=reg137/reg148; reg133=reg74*reg133; reg132=reg72*reg132; reg122=reg22*reg122;
   reg121=reg11*reg121; reg140=pow(reg140,0.5); reg136=reg136/reg141; reg128=reg128/reg141; reg135=reg66*reg135;
   reg33=reg70*reg33; reg84=reg84*reg148; reg133=reg132+reg133; reg97=reg142+reg97; reg148=reg31*reg148;
   reg11=reg43*reg128; reg22=reg41*reg136; reg80=reg31*reg80; reg122=reg121+reg122; reg124=reg43*reg124;
   reg119=reg117+reg119; reg141=reg139/reg141; reg100=reg81*reg100; reg128=reg109*reg128; reg55=reg55*reg138;
   reg37=reg36+reg37; reg135=reg33+reg135; reg138=reg81*reg138; reg144=reg144/reg140; reg125=reg41*reg125;
   reg143=reg143/reg140; reg136=reg111*reg136; reg120=reg120*reg143; reg116=reg116*reg141; reg102=reg102*reg144;
   reg140=reg145/reg140; reg80=reg122+reg80; reg143=reg89*reg143; reg128=reg136+reg128; reg11=reg22+reg11;
   reg141=reg45*reg141; reg144=reg86*reg144; reg55=reg37+reg55; reg138=reg135+reg138; reg106=reg89*reg106;
   reg100=reg119+reg100; reg130=reg86*reg130; reg124=reg125+reg124; reg148=reg97+reg148; reg104=reg45*reg104;
   reg84=reg133+reg84; reg102=reg120+reg102; reg138=reg129*reg138; reg84=reg80*reg84; reg127=reg127*reg140;
   reg55=reg100*reg55; reg130=reg106+reg130; reg115=reg65*reg115; reg144=reg143+reg144; reg140=reg65*reg140;
   reg116=reg128+reg116; reg148=reg126*reg148; reg141=reg11+reg141; reg104=reg124+reg104; reg140=reg144+reg140;
   reg127=reg102+reg127; reg138=reg55-reg138; reg148=reg84-reg148; reg115=reg130+reg115; reg141=reg131*reg141;
   reg116=reg104*reg116; reg11=0.024056261216234395431*reg148; reg22=0.13144585576580215187*reg138; reg31=0.024056261216234395431*reg138; reg33=0.024056261216234409915*reg138;
   reg36=0.13144585576580215187*reg148; reg37=0.04166666666666666908*reg148; reg41=0.035220810900864524453*reg148; reg43=0.024056261216234409915*reg148; reg45=0.035220810900864524453*reg138;
   reg140=reg134*reg140; reg127=reg115*reg127; reg141=reg116-reg141; reg55=0.04166666666666666908*reg138; reg65=reg41+reg45;
   reg66=0.024056261216234409915*reg141; reg31=reg31-reg37; reg70=0.13144585576580215187*reg141; reg45=reg45+reg36; reg11=reg11-reg55;
   reg36=reg36+reg22; reg72=0.035220810900864524453*reg141; reg33=reg37+reg33; reg37=0.024056261216234395431*reg141; reg22=reg41+reg22;
   reg140=reg127-reg140; reg41=0.04166666666666666908*reg141; reg55=reg43+reg55; reg43=0.035220810900864524453*reg140; reg45=reg70+reg45;
   reg74=0.13144585576580215187*reg140; reg36=reg36+reg72; reg22=reg72+reg22; reg70=reg65+reg70; reg66=reg31-reg66;
   reg31=0.024056261216234409915*reg140; reg11=reg11-reg41; reg65=0.04166666666666666908*reg140; reg33=reg37-reg33; reg37=0.024056261216234395431*reg140;
   reg41=reg55+reg41; Ne(0,9)+=reg66-reg65; Ne(1,10)+=reg66-reg65; Ne(2,11)+=reg66-reg65; Ne(0,12)+=reg70+reg74;
   Ne(1,13)+=reg70+reg74; Ne(2,14)+=reg70+reg74; Ne(0,6)+=reg11-reg31; Ne(1,7)+=reg11-reg31; Ne(2,8)+=reg11-reg31;
   Ne(0,15)+=reg45+reg43; Ne(1,16)+=reg45+reg43; Ne(2,17)+=reg45+reg43; Ne(0,3)+=reg33-reg65; Ne(1,4)+=reg33-reg65;
   Ne(2,5)+=reg33-reg65; Ne(0,0)+=reg37-reg41; Ne(1,1)+=reg37-reg41; Ne(2,2)+=reg37-reg41; Ne(0,18)+=reg43+reg36;
   Ne(1,19)+=reg43+reg36; Ne(2,20)+=reg43+reg36; Ne(0,21)+=reg74+reg22; Ne(1,22)+=reg74+reg22; Ne(2,23)+=reg74+reg22;

}

};
} // namespace LMT

