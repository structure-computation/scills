
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
   
   double reg0=0.21132486540518713447*e.pos(0)[0]; double reg1=0.21132486540518713447*e.pos(1)[0]; double reg2=0.21132486540518713447*e.pos(1)[1]; double reg3=0.21132486540518713447*e.pos(0)[1]; double reg4=reg2-reg3;
   double reg5=0.21132486540518713447*e.pos(1)[2]; double reg6=0.78867513459481286553*e.pos(2)[1]; double reg7=0.21132486540518713447*e.pos(0)[2]; double reg8=0.78867513459481286553*e.pos(0)[0]; double reg9=0.78867513459481286553*e.pos(0)[1];
   double reg10=reg1-reg0; double reg11=0.78867513459481286553*e.pos(2)[0]; double reg12=0.78867513459481286553*e.pos(1)[1]; double reg13=0.78867513459481286553*e.pos(1)[0]; double reg14=0.21132486540518713447*e.pos(2)[1];
   double reg15=reg5-reg7; double reg16=0.21132486540518713447*e.pos(2)[0]; double reg17=0.78867513459481286553*e.pos(1)[2]; double reg18=0.78867513459481286553*e.pos(2)[2]; double reg19=0.78867513459481286553*e.pos(0)[2];
   double reg20=0.78867513459481286553*e.pos(3)[1]; reg4=reg6+reg4; double reg21=0.78867513459481286553*e.pos(3)[0]; reg10=reg11+reg10; double reg22=reg12-reg9;
   double reg23=reg13-reg8; double reg24=0.21132486540518713447*e.pos(2)[2]; double reg25=reg17-reg19; reg22=reg22+reg14; reg23=reg16+reg23;
   double reg26=0.21132486540518713447*e.pos(3)[1]; double reg27=0.21132486540518713447*e.pos(3)[0]; reg10=reg10-reg21; reg4=reg4-reg20; reg15=reg18+reg15;
   double reg28=0.78867513459481286553*e.pos(3)[2]; double reg29=pow(reg4,2); reg23=reg23-reg27; reg22=reg22-reg26; reg15=reg15-reg28;
   double reg30=pow(reg10,2); reg25=reg24+reg25; double reg31=0.21132486540518713447*e.pos(3)[2]; double reg32=pow(reg15,2); reg25=reg25-reg31;
   reg29=reg30+reg29; reg30=pow(reg22,2); double reg33=pow(reg23,2); reg9=reg2+reg9; reg2=pow(reg25,2);
   reg32=reg29+reg32; reg8=reg1+reg8; reg12=reg3+reg12; reg30=reg33+reg30; reg13=reg0+reg13;
   reg8=reg16-reg8; reg9=reg14-reg9; reg19=reg5+reg19; reg17=reg7+reg17; reg2=reg30+reg2;
   reg12=reg6-reg12; reg32=pow(reg32,0.5); reg13=reg11-reg13; reg9=reg20+reg9; reg19=reg24-reg19;
   reg8=reg21+reg8; reg2=pow(reg2,0.5); reg0=reg4/reg32; reg1=reg10/reg32; reg27=reg13+reg27;
   reg26=reg12+reg26; reg17=reg18-reg17; reg32=reg15/reg32; reg3=reg1*reg27; reg5=reg0*reg26;
   reg19=reg28+reg19; reg31=reg17+reg31; reg6=reg0*reg9; reg7=reg1*reg8; reg11=reg23/reg2;
   reg12=reg22/reg2; reg5=reg3+reg5; reg3=reg32*reg19; reg13=reg32*reg31; reg6=reg7+reg6;
   reg7=reg27*reg11; reg2=reg25/reg2; reg14=reg26*reg12; reg16=reg9*reg12; reg17=reg8*reg11;
   reg3=reg6+reg3; reg6=reg31*reg2; reg13=reg5+reg13; reg14=reg7+reg14; reg5=reg0*reg3;
   reg7=reg1*reg3; reg6=reg14+reg6; reg16=reg17+reg16; reg14=reg19*reg2; reg17=reg0*reg13;
   reg18=reg1*reg13; reg18=reg27-reg18; reg20=reg12*reg6; reg21=reg32*reg3; reg5=reg9-reg5;
   reg24=reg32*reg13; reg7=reg8-reg7; reg14=reg16+reg14; reg17=reg26-reg17; reg16=reg11*reg6;
   reg28=pow(reg18,2); reg29=reg12*reg14; reg30=reg11*reg14; reg33=reg2*reg6; reg20=reg26-reg20;
   reg16=reg27-reg16; reg21=reg19-reg21; reg24=reg31-reg24; double reg34=pow(reg5,2); double reg35=pow(reg17,2);
   double reg36=pow(reg7,2); double reg37=pow(reg20,2); double reg38=pow(reg16,2); reg33=reg31-reg33; double reg39=reg2*reg14;
   reg30=reg8-reg30; reg29=reg9-reg29; double reg40=pow(reg24,2); reg35=reg28+reg35; reg34=reg36+reg34;
   reg28=pow(reg21,2); reg36=pow(reg29,2); reg40=reg35+reg40; reg39=reg19-reg39; reg35=pow(reg33,2);
   reg28=reg34+reg28; reg37=reg38+reg37; reg34=pow(reg30,2); reg36=reg34+reg36; reg34=pow(reg39,2);
   reg28=pow(reg28,0.5); reg35=reg37+reg35; reg40=pow(reg40,0.5); reg5=reg5/reg28; reg7=reg7/reg28;
   reg17=reg17/reg40; reg35=pow(reg35,0.5); reg34=reg36+reg34; reg18=reg18/reg40; reg40=reg24/reg40;
   reg24=reg26*reg17; reg36=reg10*reg18; reg18=reg27*reg18; reg16=reg16/reg35; reg17=reg4*reg17;
   reg20=reg20/reg35; reg0=reg4*reg0; reg37=reg10*reg7; reg4=reg4*reg5; reg34=pow(reg34,0.5);
   reg28=reg21/reg28; reg5=reg9*reg5; reg7=reg8*reg7; reg1=reg10*reg1; reg12=reg22*reg12;
   reg10=reg22*reg20; reg21=reg23*reg16; reg4=reg37+reg4; reg37=reg15*reg28; reg35=reg33/reg35;
   reg20=reg26*reg20; reg16=reg27*reg16; reg29=reg29/reg34; reg30=reg30/reg34; reg0=reg1+reg0;
   reg28=reg19*reg28; reg1=reg15*reg40; reg24=reg18+reg24; reg5=reg7+reg5; reg40=reg31*reg40;
   reg32=reg15*reg32; reg17=reg36+reg17; reg11=reg23*reg11; reg8=reg8*reg30; reg9=reg9*reg29;
   reg34=reg39/reg34; reg30=reg23*reg30; reg32=reg0+reg32; reg29=reg22*reg29; reg28=reg5+reg28;
   reg0=reg25*reg35; reg10=reg21+reg10; reg37=reg4+reg37; reg12=reg11+reg12; reg35=reg31*reg35;
   reg20=reg16+reg20; reg40=reg24+reg40; reg1=reg17+reg1; reg2=reg25*reg2; reg9=reg8+reg9;
   reg28=reg32*reg28; reg19=reg19*reg34; reg37=reg3*reg37; reg1=reg13*reg1; reg2=reg12+reg2;
   reg40=reg32*reg40; reg35=reg20+reg35; reg29=reg30+reg29; reg34=reg25*reg34; reg0=reg10+reg0;
   reg19=reg9+reg19; reg34=reg29+reg34; reg37=reg28-reg37; reg1=reg40-reg1; reg35=reg2*reg35;
   reg0=reg6*reg0; reg3=0.1555021169820365473*reg37; reg4=0.1555021169820365473*reg1; reg5=0.011164549684630114537*reg37; reg6=0.04166666666666666908*reg1;
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
   
   double reg0=2.1547005383792514621*e.pos(0)[1]; double reg1=0.15470053837925146212*e.pos(1)[1]; double reg2=2.1547005383792514621*e.pos(1)[0]; double reg3=0.15470053837925146212*e.pos(0)[0]; double reg4=2.1547005383792514621*e.pos(0)[0];
   double reg5=0.15470053837925146212*e.pos(1)[0]; double reg6=2.1547005383792514621*e.pos(1)[1]; double reg7=0.15470053837925146212*e.pos(0)[1]; double reg8=2.3094010767585029242*e.pos(2)[1]; reg6=reg6+reg7;
   reg0=reg1+reg0; reg4=reg5+reg4; double reg9=2.3094010767585029242*e.pos(2)[0]; reg2=reg2+reg3; reg6=reg6-reg8;
   reg2=reg2-reg9; double reg10=reg8-reg0; double reg11=reg9-reg4; double reg12=pow(reg6,2); double reg13=pow(reg10,2);
   double reg14=pow(reg2,2); double reg15=pow(reg11,2); reg12=reg14+reg12; reg13=reg15+reg13; reg12=pow(reg12,0.5);
   reg13=pow(reg13,0.5); reg14=reg10/reg13; reg13=reg11/reg13; reg15=reg2/reg12; reg12=reg6/reg12;
   reg14=reg10*reg14; reg13=reg11*reg13; reg15=reg2*reg15; reg12=reg6*reg12; reg12=reg15+reg12;
   reg14=reg13+reg14; reg2=0.22767090063073975644*reg14; reg6=0.061004233964073109089*reg12; reg10=0.061004233964073109089*reg14; reg11=0.22767090063073975644*reg12;
   reg13=0.33333333333333335264*reg14; reg15=0.33333333333333335264*reg12; Ne(0,0)+=reg2-reg6; Ne(1,1)+=reg2-reg6; Ne(0,2)+=reg11-reg10;
   Ne(1,3)+=reg11-reg10; Ne(0,4)+=reg13+reg15; Ne(1,5)+=reg13+reg15;

}
 /// pour les Triangle_6
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Triangle_6,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.78379396366385999995*e.pos(1)[1]; double reg1=0.56758792732771999991*e.pos(0)[1]; double reg2=0.78379396366385999995*e.pos(1)[0]; double reg3=0.56758792732771999991*e.pos(0)[0]; double reg4=0.78379396366385999995*e.pos(0)[0];
   double reg5=0.56758792732771999991*e.pos(1)[0]; double reg6=0.56758792732771999991*e.pos(1)[1]; double reg7=0.78379396366385999995*e.pos(0)[1]; double reg8=1.3513818909915799999*e.pos(3)[0]; double reg9=reg5+reg4;
   double reg10=reg0+reg1; double reg11=1.3513818909915799999*e.pos(3)[1]; double reg12=reg6+reg7; double reg13=0.78379396366385999995*e.pos(1)[2]; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg2+reg3; double reg16=0.56758792732771999991*e.pos(1)[2]; double reg17=0.78379396366385999995*e.pos(0)[2]; double reg18=1.3513818909915799999*e.pos(3)[2]; double reg19=reg13+reg14;
   double reg20=2.2673902919218320001*e.pos(0)[1]; double reg21=reg11-reg12; double reg22=2.2673902919218320001*e.pos(0)[0]; double reg23=0.63369514596091600003*e.pos(1)[0]; double reg24=1.78379396366386*e.pos(4)[1];
   reg10=reg10-reg11; double reg25=0.63369514596091600003*e.pos(1)[1]; double reg26=reg16+reg17; double reg27=reg8-reg9; double reg28=1.78379396366386*e.pos(4)[0];
   reg15=reg15-reg8; double reg29=2.2673902919218320001*e.pos(1)[0]; reg0=reg0-reg7; double reg30=0.63369514596091600003*e.pos(0)[1]; double reg31=2.2673902919218320001*e.pos(1)[1];
   double reg32=reg23+reg22; double reg33=0.63369514596091600003*e.pos(0)[0]; double reg34=2.9010854378827480001*e.pos(3)[0]; double reg35=2.2673902919218320001*e.pos(0)[2]; double reg36=0.63369514596091600003*e.pos(1)[2];
   double reg37=reg25+reg20; double reg38=2.9010854378827480001*e.pos(3)[1]; double reg39=0.43241207267228000009*e.pos(4)[0]; double reg40=0.43241207267228000009*e.pos(4)[1]; reg15=reg15+reg28;
   double reg41=1.78379396366386*e.pos(5)[0]; reg27=reg28+reg27; double reg42=reg18-reg26; double reg43=1.78379396366386*e.pos(4)[2]; reg19=reg19-reg18;
   reg21=reg24+reg21; double reg44=1.78379396366386*e.pos(5)[1]; reg10=reg10+reg24; reg2=reg2-reg4; double reg45=0.366304854039084*e.pos(4)[0];
   reg42=reg43+reg42; double reg46=0.43241207267228000009*e.pos(4)[2]; double reg47=reg38-reg37; double reg48=reg34-reg32; reg21=reg21-reg44;
   reg27=reg27-reg41; reg2=reg39+reg2; reg13=reg13-reg17; double reg49=0.43241207267228000009*e.pos(5)[1]; reg0=reg40+reg0;
   double reg50=0.43241207267228000009*e.pos(5)[0]; reg41=reg15-reg41; reg44=reg10-reg44; reg19=reg19+reg43; reg10=1.78379396366386*e.pos(5)[2];
   reg15=0.63369514596091600003*e.pos(0)[2]; double reg51=2.2673902919218320001*e.pos(1)[2]; reg31=reg31+reg30; reg29=reg29+reg33; double reg52=0.366304854039084*e.pos(4)[1];
   double reg53=2.9010854378827480001*e.pos(3)[2]; double reg54=reg36+reg35; double reg55=reg53-reg54; reg48=reg48+reg45; reg51=reg51+reg15;
   double reg56=2.710505431213761085e-20*e.pos(3)[1]; reg50=reg2-reg50; reg2=pow(reg21,2); double reg57=pow(reg41,2); double reg58=reg33-reg23;
   double reg59=0.366304854039084*e.pos(5)[0]; double reg60=2.710505431213761085e-20*e.pos(3)[0]; reg19=reg19-reg10; double reg61=pow(reg27,2); double reg62=pow(reg44,2);
   double reg63=reg30-reg25; double reg64=0.78379396366385999995*e.pos(2)[1]; reg29=reg29-reg34; reg47=reg47+reg52; double reg65=0.366304854039084*e.pos(4)[2];
   double reg66=0.78379396366385999995*e.pos(2)[0]; reg31=reg31-reg38; reg10=reg42-reg10; reg42=0.366304854039084*e.pos(5)[1]; double reg67=0.43241207267228000009*e.pos(5)[2];
   reg13=reg46+reg13; reg49=reg0-reg49; reg47=reg47-reg42; reg0=reg15-reg36; double reg68=2.710505431213761085e-20*e.pos(3)[2];
   reg63=reg63-reg56; reg48=reg48-reg59; reg58=reg58-reg60; double reg69=3.2673902919218320001*e.pos(4)[1]; double reg70=3.2673902919218320001*e.pos(4)[0];
   reg51=reg51-reg53; double reg71=pow(reg50,2); double reg72=pow(reg49,2); reg31=reg52+reg31; reg67=reg13-reg67;
   reg13=0.56758792732771999991*e.pos(2)[0]; double reg73=0.56758792732771999991*e.pos(2)[1]; reg29=reg45+reg29; reg55=reg55+reg65; double reg74=0.366304854039084*e.pos(5)[2];
   double reg75=pow(reg19,2); double reg76=pow(reg10,2); double reg77=reg3+reg66; double reg78=1.78379396366386*e.pos(3)[0]; double reg79=1.78379396366386*e.pos(3)[1];
   double reg80=reg1+reg64; double reg81=0.78379396366385999995*e.pos(2)[2]; reg2=reg61+reg2; reg62=reg57+reg62; reg58=reg70+reg58;
   reg57=3.2673902919218320001*e.pos(5)[0]; reg63=reg69+reg63; reg61=3.2673902919218320001*e.pos(5)[1]; double reg82=0.63369514596091600003*e.pos(2)[1]; reg0=reg0-reg68;
   double reg83=0.63369514596091600003*e.pos(2)[0]; reg75=reg62+reg75; reg62=3.2673902919218320001*e.pos(4)[2]; reg59=reg29-reg59; reg29=0.56758792732771999991*e.pos(2)[2];
   reg77=reg77-reg78; reg72=reg71+reg72; reg71=reg7+reg73; double reg84=pow(reg67,2); double reg85=reg4+reg13;
   reg51=reg65+reg51; reg42=reg31-reg42; reg31=1.78379396366386*e.pos(3)[2]; double reg86=reg14+reg81; double reg87=pow(reg47,2);
   reg7=reg64-reg7; reg55=reg55-reg74; reg64=pow(reg48,2); reg4=reg66-reg4; reg66=0.43241207267228000009*e.pos(3)[1];
   reg80=reg80-reg79; double reg88=0.43241207267228000009*e.pos(3)[0]; reg76=reg2+reg76; reg2=0.63369514596091600003*e.pos(2)[2]; reg74=reg51-reg74;
   reg76=pow(reg76,0.5); reg78=reg78+reg85; reg77=reg28+reg77; reg88=reg4-reg88; reg4=reg17+reg29;
   reg51=pow(reg59,2); reg84=reg72+reg84; reg87=reg64+reg87; reg79=reg79+reg71; reg86=reg86-reg31;
   reg64=pow(reg42,2); reg22=reg22+reg83; reg72=0.366304854039084*e.pos(3)[0]; double reg89=3.2673902919218320001*e.pos(5)[2]; reg0=reg62+reg0;
   double reg90=pow(reg55,2); double reg91=0.43241207267228000009*e.pos(3)[2]; reg61=reg63-reg61; reg17=reg81-reg17; reg80=reg24+reg80;
   reg20=reg20+reg82; reg57=reg58-reg57; reg58=1.3513818909915799999*e.pos(5)[1]; reg66=reg7-reg66; reg75=pow(reg75,0.5);
   reg7=0.366304854039084*e.pos(3)[1]; reg63=1.3513818909915799999*e.pos(5)[0]; reg28=reg28-reg78; reg86=reg43+reg86; reg81=1.3513818909915799999*e.pos(5)[2];
   reg24=reg24-reg79; reg31=reg31+reg4; reg77=reg77-reg63; reg80=reg80-reg58; double reg92=0.366304854039084*e.pos(3)[2];
   reg35=reg35+reg2; double reg93=reg20+reg7; double reg94=reg22+reg72; reg90=reg87+reg90; reg87=reg27/reg76;
   double reg95=reg21/reg76; double reg96=2.2673902919218320001*e.pos(2)[1]; reg39=reg88+reg39; reg40=reg66+reg40; reg91=reg17-reg91;
   reg17=2.2673902919218320001*e.pos(2)[0]; reg89=reg0-reg89; reg0=pow(reg61,2); reg66=pow(reg57,2); reg88=reg41/reg75;
   double reg97=3.2673902919218320001*e.pos(3)[1]; double reg98=reg30-reg82; double reg99=reg44/reg75; reg64=reg51+reg64; reg51=3.2673902919218320001*e.pos(3)[0];
   double reg100=reg33-reg83; double reg101=pow(reg74,2); reg84=pow(reg84,0.5); reg0=reg66+reg0; reg101=reg64+reg101;
   reg64=pow(reg89,2); reg66=2.9010854378827480001*e.pos(5)[0]; reg90=pow(reg90,0.5); reg17=reg33+reg17; reg46=reg91+reg46;
   reg28=reg63+reg28; reg91=reg95*reg40; double reg102=reg87*reg39; double reg103=reg99*reg80; double reg104=reg50/reg84;
   double reg105=reg49/reg84; reg96=reg30+reg96; reg76=reg10/reg76; double reg106=2.2673902919218320001*e.pos(2)[2]; reg75=reg19/reg75;
   reg86=reg86-reg81; reg51=reg100-reg51; reg43=reg43-reg31; reg24=reg58+reg24; reg97=reg98-reg97;
   reg98=reg88*reg77; reg100=reg35+reg92; double reg107=reg15-reg2; double reg108=3.2673902919218320001*e.pos(3)[2]; double reg109=2.9010854378827480001*e.pos(5)[1];
   double reg110=reg45-reg94; double reg111=reg52-reg93; double reg112=reg75*reg86; reg70=reg51+reg70; reg106=reg15+reg106;
   reg110=reg66+reg110; reg101=pow(reg101,0.5); reg7=reg96-reg7; reg108=reg107-reg108; reg64=reg0+reg64;
   reg69=reg97+reg69; reg0=2.9010854378827480001*e.pos(5)[2]; reg103=reg98+reg103; reg51=reg65-reg100; reg96=5.42101086242752217e-20*e.pos(5)[0];
   reg97=5.42101086242752217e-20*e.pos(5)[1]; reg98=reg47/reg90; reg72=reg17-reg72; reg17=reg48/reg90; reg111=reg111+reg109;
   reg91=reg102+reg91; reg102=reg76*reg46; reg84=reg67/reg84; reg107=reg104*reg28; double reg113=reg105*reg24;
   reg43=reg81+reg43; reg72=reg45+reg72; reg7=reg52+reg7; reg45=reg17*reg110; reg102=reg91+reg102;
   reg90=reg55/reg90; reg52=reg98*reg111; reg92=reg106-reg92; reg51=reg51+reg0; reg64=pow(reg64,0.5);
   reg91=reg84*reg43; reg70=reg70-reg96; reg112=reg103+reg112; reg69=reg69-reg97; reg62=reg108+reg62;
   reg103=5.42101086242752217e-20*e.pos(5)[2]; reg106=reg42/reg101; reg108=reg59/reg101; reg113=reg107+reg113; reg107=reg61/reg64;
   double reg114=reg90*reg51; double reg115=reg57/reg64; reg52=reg45+reg52; reg101=reg74/reg101; reg62=reg62-reg103;
   reg45=reg87*reg102; double reg116=reg106*reg69; double reg117=reg95*reg102; double reg118=reg108*reg70; reg91=reg113+reg91;
   reg113=reg88*reg112; double reg119=reg99*reg112; reg92=reg65+reg92; reg72=reg72-reg66; reg7=reg7-reg109;
   reg113=reg77-reg113; reg65=reg115*reg72; reg116=reg118+reg116; reg119=reg80-reg119; reg118=reg101*reg62;
   double reg120=reg75*reg112; double reg121=reg104*reg91; reg117=reg40-reg117; double reg122=reg76*reg102; reg45=reg39-reg45;
   reg64=reg89/reg64; reg114=reg52+reg114; reg92=reg92-reg0; reg52=reg107*reg7; double reg123=reg105*reg91;
   double reg124=pow(reg45,2); double reg125=reg84*reg91; reg122=reg46-reg122; reg121=reg28-reg121; double reg126=pow(reg117,2);
   reg120=reg86-reg120; reg118=reg116+reg118; reg116=pow(reg119,2); double reg127=reg98*reg114; double reg128=pow(reg113,2);
   reg52=reg65+reg52; reg65=reg17*reg114; reg123=reg24-reg123; double reg129=reg64*reg92; reg116=reg128+reg116;
   reg128=reg108*reg118; double reg130=pow(reg120,2); double reg131=reg106*reg118; reg125=reg43-reg125; double reg132=pow(reg122,2);
   reg124=reg126+reg124; reg126=pow(reg121,2); reg65=reg110-reg65; reg127=reg111-reg127; reg129=reg52+reg129;
   reg52=reg90*reg114; double reg133=pow(reg123,2); double reg134=pow(reg65,2); double reg135=pow(reg127,2); double reg136=pow(reg125,2);
   reg52=reg51-reg52; reg128=reg70-reg128; reg131=reg69-reg131; double reg137=reg101*reg118; reg130=reg116+reg130;
   reg133=reg126+reg133; reg116=reg115*reg129; reg126=reg107*reg129; reg132=reg124+reg132; reg116=reg72-reg116;
   reg124=reg64*reg129; reg136=reg133+reg136; reg135=reg134+reg135; reg132=pow(reg132,0.5); reg133=pow(reg52,2);
   reg134=pow(reg128,2); double reg138=pow(reg131,2); reg137=reg62-reg137; reg126=reg7-reg126; reg130=pow(reg130,0.5);
   double reg139=pow(reg116,2); double reg140=pow(reg126,2); reg119=reg119/reg130; reg136=pow(reg136,0.5); reg113=reg113/reg130;
   reg124=reg92-reg124; reg45=reg45/reg132; reg133=reg135+reg133; reg117=reg117/reg132; reg138=reg134+reg138;
   reg134=pow(reg137,2); reg39=reg39*reg45; reg45=reg27*reg45; reg132=reg122/reg132; reg40=reg40*reg117;
   reg117=reg21*reg117; reg140=reg139+reg140; reg122=pow(reg124,2); reg95=reg21*reg95; reg87=reg27*reg87;
   reg88=reg41*reg88; reg134=reg138+reg134; reg99=reg44*reg99; reg130=reg120/reg130; reg41=reg41*reg113;
   reg44=reg44*reg119; reg133=pow(reg133,0.5); reg121=reg121/reg136; reg123=reg123/reg136; reg113=reg77*reg113;
   reg119=reg80*reg119; reg99=reg88+reg99; reg75=reg19*reg75; reg134=pow(reg134,0.5); reg117=reg45+reg117;
   reg122=reg140+reg122; reg46=reg46*reg132; reg40=reg39+reg40; reg119=reg113+reg119; reg21=reg49*reg123;
   reg27=reg50*reg121; reg127=reg127/reg133; reg65=reg65/reg133; reg95=reg87+reg95; reg76=reg10*reg76;
   reg121=reg28*reg121; reg136=reg125/reg136; reg123=reg24*reg123; reg132=reg10*reg132; reg44=reg41+reg44;
   reg19=reg19*reg130; reg130=reg86*reg130; reg105=reg49*reg105; reg104=reg50*reg104; reg84=reg67*reg84;
   reg46=reg40+reg46; reg75=reg99+reg75; reg10=reg47*reg127; reg24=reg48*reg65; reg133=reg52/reg133;
   reg127=reg111*reg127; reg65=reg110*reg65; reg132=reg117+reg132; reg67=reg67*reg136; reg19=reg44+reg19;
   reg123=reg121+reg123; reg136=reg43*reg136; reg76=reg95+reg76; reg98=reg47*reg98; reg17=reg48*reg17;
   reg21=reg27+reg21; reg105=reg104+reg105; reg122=pow(reg122,0.5); reg130=reg119+reg130; reg128=reg128/reg134;
   reg131=reg131/reg134; reg126=reg126/reg122; reg51=reg51*reg133; reg127=reg65+reg127; reg27=reg42*reg131;
   reg28=reg59*reg128; reg132=reg102*reg132; reg134=reg137/reg134; reg131=reg69*reg131; reg19=reg112*reg19;
   reg136=reg123+reg136; reg90=reg55*reg90; reg98=reg17+reg98; reg128=reg70*reg128; reg67=reg21+reg67;
   reg106=reg42*reg106; reg116=reg116/reg122; reg46=reg76*reg46; reg84=reg105+reg84; reg108=reg59*reg108;
   reg130=reg75*reg130; reg133=reg55*reg133; reg10=reg24+reg10; reg67=reg91*reg67; reg19=reg130-reg19;
   reg72=reg72*reg116; reg90=reg98+reg90; reg136=reg84*reg136; reg106=reg108+reg106; reg101=reg74*reg101;
   reg116=reg57*reg116; reg7=reg7*reg126; reg51=reg127+reg51; reg74=reg74*reg134; reg133=reg10+reg133;
   reg27=reg28+reg27; reg115=reg57*reg115; reg107=reg61*reg107; reg132=reg46-reg132; reg122=reg124/reg122;
   reg134=reg62*reg134; reg126=reg61*reg126; reg131=reg128+reg131; reg10=0.088847818743090689935*reg19; reg126=reg116+reg126;
   reg7=reg72+reg7; reg17=reg89*reg122; reg21=0.005384432036113586778*reg19; reg24=0.088847818743090689935*reg132; reg122=reg92*reg122;
   reg28=0.021537728144454347112*reg132; reg39=0.021537728144454347112*reg19; reg40=0.009463616120767210603*reg132; reg64=reg89*reg64; reg107=reg115+reg107;
   reg41=0.005384432036113586778*reg132; reg42=0.009463616120767210603*reg19; reg133=reg114*reg133; reg67=reg136-reg67; reg51=reg90*reg51;
   reg74=reg27+reg74; reg134=reg131+reg134; reg101=reg106+reg101; reg24=reg39+reg24; reg40=reg21+reg40;
   reg21=reg41+reg21; reg27=0.009463616120767210603*reg67; reg43=0.021537728144454347112*reg67; reg10=reg28+reg10; reg41=reg42+reg41;
   reg42=0.005384432036113586778*reg67; reg28=reg39+reg28; reg39=0.088847818743090689935*reg67; reg64=reg107+reg64; reg17=reg126+reg17;
   reg133=reg51-reg133; reg74=reg118*reg74; reg134=reg101*reg134; reg122=reg7+reg122; reg7=0.016449618187943419918*reg133;
   reg39=reg28+reg39; reg74=reg134-reg74; reg41=reg41+reg42; reg122=reg64*reg122; reg27=reg21+reg27;
   reg10=reg10+reg43; reg21=0.028457289286966203713*reg133; reg28=0.0018441552587796664112*reg133; reg44=0.0041124045469858549794*reg133; reg40=reg42+reg40;
   reg17=reg129*reg17; reg24=reg43+reg24; reg24=reg7+reg24; reg28=reg10+reg28; reg10=0.0018441552587796664109*reg74;
   reg42=0.016449618187943419918*reg74; reg43=0.016449618187943419916*reg74; reg45=0.004112404546985854979*reg74; reg41=reg21-reg41; reg17=reg122-reg17;
   reg21=0.028457289286966203713*reg74; reg40=reg40+reg44; reg27=reg44+reg27; reg44=0.0041124045469858549794*reg74; reg7=reg39+reg7;
   reg39=0.028457289286966203713*reg17; reg44=reg27+reg44; reg40=reg21-reg40; reg10=reg24+reg10; reg21=0.0018441552587796664111*reg17;
   reg42=reg28+reg42; reg24=0.016449618187943419918*reg17; reg45=reg41-reg45; reg27=0.0041124045469858549794*reg17; reg43=reg7+reg43;
   Ne(0,15)+=reg24+reg10; Ne(1,16)+=reg24+reg10; Ne(2,17)+=reg24+reg10; Ne(0,6)+=reg39-reg44; Ne(1,7)+=reg39-reg44;
   Ne(2,8)+=reg39-reg44; Ne(0,0)+=reg45-reg27; Ne(1,1)+=reg45-reg27; Ne(2,2)+=reg45-reg27; Ne(0,12)+=reg42+reg24;
   Ne(1,13)+=reg42+reg24; Ne(2,14)+=reg42+reg24; Ne(0,9)+=reg43+reg21; Ne(1,10)+=reg43+reg21; Ne(2,11)+=reg43+reg21;
   Ne(0,3)+=reg40-reg27; Ne(1,4)+=reg40-reg27; Ne(2,5)+=reg40-reg27;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(0)[0]; double reg1=0.12200846792814621817*e.pos(0)[0]; double reg2=0.36602540378443865451*e.pos(1)[0]; double reg3=0.12200846792814621817*e.pos(1)[1]; double reg4=0.36602540378443865451*e.pos(0)[1];
   double reg5=0.12200846792814621817*e.pos(1)[0]; double reg6=0.36602540378443865451*e.pos(1)[1]; double reg7=0.12200846792814621817*e.pos(0)[1]; double reg8=0.45534180126147951289*e.pos(0)[0]; double reg9=0.45534180126147951289*e.pos(0)[1];
   double reg10=0.12200846792814621817*e.pos(0)[2]; double reg11=0.36602540378443865451*e.pos(1)[2]; double reg12=reg5+reg0; double reg13=1.3660254037844385386*e.pos(1)[0]; double reg14=1.3660254037844385386*e.pos(1)[1];
   double reg15=1.3660254037844385386*e.pos(2)[0]; double reg16=0.45534180126147951289*e.pos(2)[1]; double reg17=reg6+reg7; double reg18=0.36602540378443865451*e.pos(0)[2]; double reg19=0.12200846792814621817*e.pos(1)[2];
   double reg20=1.3660254037844385386*e.pos(2)[1]; double reg21=reg3+reg4; double reg22=0.45534180126147951289*e.pos(2)[0]; double reg23=reg2+reg1; double reg24=reg19+reg18;
   double reg25=1.3660254037844385386*e.pos(2)[2]; double reg26=1.3660254037844385386*e.pos(1)[2]; double reg27=reg17+reg16; double reg28=reg9+reg14; double reg29=0.45534180126147951289*e.pos(1)[0];
   double reg30=1.3660254037844385386*e.pos(3)[1]; double reg31=reg8+reg13; double reg32=0.45534180126147951289*e.pos(1)[1]; double reg33=1.3660254037844385386*e.pos(0)[1]; double reg34=reg11+reg10;
   double reg35=0.45534180126147951289*e.pos(2)[2]; double reg36=0.45534180126147951289*e.pos(0)[2]; double reg37=0.12200846792814621817*e.pos(2)[1]; double reg38=0.12200846792814621817*e.pos(2)[0]; double reg39=reg23+reg22;
   double reg40=0.45534180126147951289*e.pos(3)[0]; double reg41=1.3660254037844385386*e.pos(0)[0]; reg12=reg12+reg15; double reg42=1.3660254037844385386*e.pos(3)[0]; reg21=reg21+reg20;
   double reg43=0.45534180126147951289*e.pos(3)[1]; reg12=reg12+reg40; double reg44=0.36602540378443865451*e.pos(2)[0]; double reg45=0.36602540378443865451*e.pos(2)[1]; double reg46=reg36+reg26;
   double reg47=reg39+reg42; double reg48=reg32+reg33; double reg49=reg27+reg30; double reg50=0.36602540378443865451*e.pos(3)[1]; reg28=reg37+reg28;
   double reg51=0.48803387171258487271*e.pos(4)[0]; double reg52=0.48803387171258487271*e.pos(4)[1]; double reg53=reg29+reg41; double reg54=0.12200846792814621817*e.pos(2)[2]; double reg55=0.36602540378443865451*e.pos(3)[0];
   double reg56=1.3660254037844385386*e.pos(3)[2]; double reg57=reg34+reg35; reg21=reg21+reg43; reg24=reg24+reg25; double reg58=0.45534180126147951289*e.pos(3)[2];
   reg31=reg38+reg31; double reg59=0.45534180126147951289*e.pos(1)[2]; double reg60=1.3660254037844385386*e.pos(0)[2]; reg12=reg12-reg51; reg31=reg31+reg55;
   double reg61=reg53+reg44; double reg62=0.12200846792814621817*e.pos(3)[1]; double reg63=reg48+reg45; double reg64=reg59+reg60; double reg65=0.36602540378443865451*e.pos(2)[2];
   double reg66=reg57+reg56; double reg67=0.66666666666666670528*e.pos(5)[0]; reg21=reg21-reg52; double reg68=0.66666666666666670528*e.pos(5)[1]; double reg69=reg51-reg47;
   reg24=reg24+reg58; double reg70=0.48803387171258487271*e.pos(4)[2]; double reg71=0.36602540378443865451*e.pos(3)[2]; reg46=reg54+reg46; double reg72=1.8213672050459180516*e.pos(4)[1];
   double reg73=1.8213672050459180516*e.pos(4)[0]; double reg74=reg52-reg49; double reg75=0.12200846792814621817*e.pos(3)[0]; reg28=reg28+reg50; reg31=reg31-reg73;
   reg21=reg21+reg68; double reg76=1.8213672050459180516*e.pos(6)[1]; double reg77=reg70-reg66; reg74=reg68+reg74; reg69=reg67+reg69;
   double reg78=1.8213672050459180516*e.pos(4)[2]; reg46=reg46+reg71; double reg79=reg64+reg65; double reg80=0.12200846792814621817*e.pos(3)[2]; double reg81=reg75+reg61;
   double reg82=1.8213672050459180516*e.pos(6)[0]; reg24=reg24-reg70; double reg83=0.66666666666666670528*e.pos(5)[2]; double reg84=reg62+reg63; reg12=reg67+reg12;
   reg28=reg28-reg72; reg69=reg82+reg69; reg77=reg83+reg77; reg74=reg76+reg74; reg12=reg12-reg82;
   double reg85=0.66666666666666670528*e.pos(7)[0]; reg21=reg21-reg76; double reg86=0.66666666666666670528*e.pos(7)[1]; reg46=reg46-reg78; reg24=reg24+reg83;
   double reg87=1.8213672050459180516*e.pos(6)[2]; double reg88=0.48803387171258487271*e.pos(6)[1]; reg28=reg68+reg28; double reg89=reg0+reg29; double reg90=0.48803387171258487271*e.pos(6)[0];
   reg31=reg67+reg31; double reg91=reg4+reg32; double reg92=reg73-reg81; double reg93=reg2+reg8; double reg94=reg72-reg84;
   double reg95=reg80+reg79; double reg96=reg6+reg9; reg14=reg7+reg14; reg13=reg1+reg13; reg38=reg38+reg93;
   reg77=reg87+reg77; reg1=0.48803387171258487271*e.pos(6)[2]; reg46=reg83+reg46; reg37=reg37+reg96; reg7=reg11+reg36;
   reg69=reg69-reg85; double reg97=reg78-reg95; reg28=reg28-reg88; reg94=reg68+reg94; reg92=reg67+reg92;
   reg31=reg31-reg90; reg74=reg74-reg86; reg67=reg18+reg59; reg12=reg12-reg85; reg20=reg20+reg91;
   reg21=reg21-reg86; reg15=reg15+reg89; reg68=0.66666666666666670528*e.pos(7)[2]; reg24=reg24-reg87; reg13=reg22+reg13;
   reg28=reg28-reg86; reg54=reg54+reg7; reg31=reg31-reg85; reg75=reg15+reg75; reg77=reg77-reg68;
   reg15=0.66666666666666670528*e.pos(4)[0]; reg38=reg42+reg38; reg62=reg20+reg62; reg20=0.66666666666666670528*e.pos(4)[1]; reg14=reg16+reg14;
   reg97=reg83+reg97; reg92=reg90+reg92; reg25=reg25+reg67; reg46=reg46-reg1; reg41=reg5+reg41;
   reg5=pow(reg69,2); reg42=pow(reg12,2); reg83=pow(reg21,2); reg33=reg3+reg33; reg37=reg30+reg37;
   reg24=reg24-reg68; reg3=pow(reg74,2); reg94=reg88+reg94; reg26=reg10+reg26; reg80=reg25+reg80;
   reg10=reg55+reg13; reg37=reg37-reg20; reg85=reg92-reg85; reg25=0.48803387171258487271*e.pos(5)[0]; reg38=reg38-reg15;
   reg30=1.8213672050459180516*e.pos(5)[1]; reg86=reg94-reg86; reg92=0.48803387171258487271*e.pos(5)[1]; reg94=0.66666666666666670528*e.pos(4)[2]; double reg98=reg50+reg14;
   reg46=reg46-reg68; reg41=reg44+reg41; reg83=reg42+reg83; reg33=reg45+reg33; reg42=pow(reg24,2);
   reg54=reg56+reg54; reg3=reg5+reg3; reg26=reg35+reg26; reg5=pow(reg28,2); reg60=reg19+reg60;
   reg19=pow(reg31,2); reg56=pow(reg77,2); reg75=reg75-reg15; double reg99=1.8213672050459180516*e.pos(5)[0]; reg97=reg1+reg97;
   reg62=reg62-reg20; double reg100=reg15+reg10; reg80=reg80-reg94; double reg101=0.48803387171258487271*e.pos(5)[2]; double reg102=pow(reg86,2);
   double reg103=0.66666666666666670528*e.pos(6)[1]; reg62=reg62-reg30; double reg104=pow(reg85,2); double reg105=0.66666666666666670528*e.pos(6)[0]; reg75=reg75-reg99;
   reg42=reg83+reg42; reg5=reg19+reg5; reg19=pow(reg46,2); reg60=reg65+reg60; reg37=reg37-reg92;
   reg38=reg38-reg25; reg83=reg20+reg98; double reg106=reg43+reg33; double reg107=reg71+reg26; reg68=reg97-reg68;
   reg97=1.8213672050459180516*e.pos(5)[2]; double reg108=reg40+reg41; reg56=reg3+reg56; reg54=reg54-reg94; reg3=pow(reg68,2);
   double reg109=reg30-reg83; double reg110=reg58+reg60; double reg111=0.48803387171258487271*e.pos(7)[0]; reg75=reg75+reg105; reg38=reg105+reg38;
   reg56=pow(reg56,0.5); double reg112=1.8213672050459180516*e.pos(7)[0]; reg62=reg62+reg103; double reg113=0.48803387171258487271*e.pos(7)[1]; reg19=reg5+reg19;
   reg54=reg54-reg101; reg15=reg15+reg108; reg5=reg94+reg107; double reg114=reg99-reg100; reg80=reg80-reg97;
   double reg115=0.66666666666666670528*e.pos(6)[2]; reg20=reg20+reg106; double reg116=1.8213672050459180516*e.pos(7)[1]; reg37=reg103+reg37; reg102=reg104+reg102;
   reg42=pow(reg42,0.5); reg104=reg92-reg20; reg3=reg102+reg3; reg94=reg94+reg110; reg19=pow(reg19,0.5);
   reg109=reg103+reg109; reg102=reg25-reg15; double reg117=reg21/reg42; double reg118=reg97-reg5; double reg119=reg12/reg42;
   reg80=reg80+reg115; double reg120=0.48803387171258487271*e.pos(7)[2]; reg37=reg37-reg116; reg62=reg62-reg113; reg114=reg105+reg114;
   reg38=reg38-reg112; double reg121=1.8213672050459180516*e.pos(7)[2]; reg54=reg115+reg54; reg75=reg75-reg111; double reg122=reg74/reg56;
   double reg123=reg69/reg56; reg54=reg54-reg121; reg114=reg111+reg114; reg102=reg105+reg102; reg105=reg122*reg37;
   reg104=reg103+reg104; reg118=reg115+reg118; reg109=reg113+reg109; reg80=reg80-reg120; reg103=reg117*reg62;
   double reg124=reg123*reg38; double reg125=reg31/reg19; reg42=reg24/reg42; double reg126=reg28/reg19; double reg127=reg101-reg94;
   reg56=reg77/reg56; double reg128=reg119*reg75; reg3=pow(reg3,0.5); double reg129=reg126*reg109; reg19=reg46/reg19;
   double reg130=reg125*reg114; reg118=reg120+reg118; double reg131=reg42*reg80; reg105=reg124+reg105; reg124=reg56*reg54;
   reg127=reg115+reg127; reg115=reg85/reg3; reg103=reg128+reg103; reg128=reg86/reg3; reg104=reg116+reg104;
   reg102=reg112+reg102; reg131=reg103+reg131; reg124=reg105+reg124; reg103=reg128*reg104; reg105=reg115*reg102;
   reg3=reg68/reg3; double reg132=reg19*reg118; reg127=reg121+reg127; reg129=reg130+reg129; reg130=reg119*reg131;
   double reg133=reg117*reg131; reg103=reg105+reg103; reg105=reg123*reg124; double reg134=reg122*reg124; double reg135=reg3*reg127;
   reg132=reg129+reg132; reg129=reg125*reg132; double reg136=reg42*reg131; reg133=reg62-reg133; reg130=reg75-reg130;
   reg105=reg38-reg105; reg134=reg37-reg134; double reg137=reg56*reg124; double reg138=reg126*reg132; reg135=reg103+reg135;
   reg138=reg109-reg138; reg103=reg19*reg132; double reg139=reg128*reg135; double reg140=reg115*reg135; reg137=reg54-reg137;
   double reg141=pow(reg134,2); double reg142=pow(reg105,2); reg136=reg80-reg136; double reg143=pow(reg133,2); double reg144=pow(reg130,2);
   reg129=reg114-reg129; reg103=reg118-reg103; double reg145=pow(reg136,2); reg143=reg144+reg143; reg144=pow(reg137,2);
   double reg146=pow(reg138,2); double reg147=reg3*reg135; reg139=reg104-reg139; double reg148=pow(reg129,2); reg140=reg102-reg140;
   reg141=reg142+reg141; reg145=reg143+reg145; reg146=reg148+reg146; reg142=pow(reg140,2); reg143=pow(reg139,2);
   reg147=reg127-reg147; reg144=reg141+reg144; reg141=pow(reg103,2); reg145=pow(reg145,0.5); reg144=pow(reg144,0.5);
   reg148=pow(reg147,2); reg141=reg146+reg141; reg143=reg142+reg143; reg134=reg134/reg144; reg105=reg105/reg144;
   reg141=pow(reg141,0.5); reg133=reg133/reg145; reg130=reg130/reg145; reg148=reg143+reg148; reg142=reg12*reg130;
   reg143=reg21*reg133; reg123=reg69*reg123; reg122=reg74*reg122; reg148=pow(reg148,0.5); reg117=reg21*reg117;
   reg119=reg12*reg119; reg145=reg136/reg145; reg133=reg62*reg133; reg130=reg75*reg130; reg38=reg38*reg105;
   reg138=reg138/reg141; reg129=reg129/reg141; reg74=reg74*reg134; reg105=reg69*reg105; reg144=reg137/reg144;
   reg134=reg37*reg134; reg42=reg24*reg42; reg117=reg119+reg117; reg143=reg142+reg143; reg24=reg24*reg145;
   reg141=reg103/reg141; reg145=reg80*reg145; reg109=reg109*reg138; reg114=reg114*reg129; reg134=reg38+reg134;
   reg129=reg31*reg129; reg133=reg130+reg133; reg140=reg140/reg148; reg125=reg31*reg125; reg139=reg139/reg148;
   reg138=reg28*reg138; reg122=reg123+reg122; reg56=reg77*reg56; reg126=reg28*reg126; reg77=reg77*reg144;
   reg74=reg105+reg74; reg144=reg54*reg144; reg42=reg117+reg42; reg138=reg129+reg138; reg12=reg46*reg141;
   reg141=reg118*reg141; reg109=reg114+reg109; reg126=reg125+reg126; reg19=reg46*reg19; reg102=reg102*reg140;
   reg104=reg104*reg139; reg148=reg147/reg148; reg140=reg85*reg140; reg139=reg86*reg139; reg77=reg74+reg77;
   reg115=reg85*reg115; reg128=reg86*reg128; reg145=reg133+reg145; reg56=reg122+reg56; reg144=reg134+reg144;
   reg24=reg143+reg24; reg77=reg124*reg77; reg19=reg126+reg19; reg104=reg102+reg104; reg128=reg115+reg128;
   reg3=reg68*reg3; reg127=reg127*reg148; reg24=reg131*reg24; reg144=reg56*reg144; reg141=reg109+reg141;
   reg139=reg140+reg139; reg145=reg42*reg145; reg148=reg68*reg148; reg12=reg138+reg12; reg127=reg104+reg127;
   reg148=reg139+reg148; reg24=reg145-reg24; reg3=reg128+reg3; reg77=reg144-reg77; reg141=reg19*reg141;
   reg12=reg132*reg12; reg19=0.024056261216234395431*reg77; reg21=0.024056261216234395431*reg24; reg28=0.035220810900864524453*reg24; reg31=0.035220810900864524453*reg77;
   reg37=0.024056261216234409915*reg77; reg12=reg141-reg12; reg38=0.04166666666666666908*reg24; reg42=0.13144585576580215187*reg24; reg148=reg135*reg148;
   reg46=0.13144585576580215187*reg77; reg54=0.04166666666666666908*reg77; reg127=reg3*reg127; reg3=0.024056261216234409915*reg24; reg56=0.13144585576580215187*reg12;
   reg62=reg42+reg46; reg42=reg31+reg42; reg31=reg28+reg31; reg68=0.035220810900864524453*reg12; reg69=0.024056261216234409915*reg12;
   reg19=reg19-reg38; reg74=0.04166666666666666908*reg12; reg3=reg3+reg54; reg148=reg127-reg148; reg75=0.024056261216234395431*reg12;
   reg37=reg38+reg37; reg46=reg28+reg46; reg54=reg21-reg54; reg3=reg3+reg74; reg21=0.024056261216234409915*reg148;
   reg28=0.035220810900864524453*reg148; reg42=reg56+reg42; reg74=reg54-reg74; reg38=0.024056261216234395431*reg148; reg62=reg62+reg68;
   reg54=0.13144585576580215187*reg148; reg56=reg31+reg56; reg69=reg19-reg69; reg37=reg75-reg37; reg19=0.04166666666666666908*reg148;
   reg46=reg68+reg46; Ne(0,21)+=reg54+reg46; Ne(1,22)+=reg54+reg46; Ne(2,23)+=reg54+reg46; Ne(0,18)+=reg28+reg62;
   Ne(1,19)+=reg28+reg62; Ne(2,20)+=reg28+reg62; Ne(0,15)+=reg42+reg28; Ne(1,16)+=reg42+reg28; Ne(2,17)+=reg42+reg28;
   Ne(0,12)+=reg56+reg54; Ne(1,13)+=reg56+reg54; Ne(2,14)+=reg56+reg54; Ne(0,9)+=reg69-reg19; Ne(1,10)+=reg69-reg19;
   Ne(2,11)+=reg69-reg19; Ne(0,6)+=reg74-reg21; Ne(1,7)+=reg74-reg21; Ne(2,8)+=reg74-reg21; Ne(0,3)+=reg37-reg19;
   Ne(1,4)+=reg37-reg19; Ne(2,5)+=reg37-reg19; Ne(0,0)+=reg38-reg3; Ne(1,1)+=reg38-reg3; Ne(2,2)+=reg38-reg3;

}

};
} // namespace LMT

