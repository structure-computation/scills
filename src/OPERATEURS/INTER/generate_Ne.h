
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
   
   double reg0=0.15470053837925146212*e.pos(1)[1]; double reg1=2.1547005383792514621*e.pos(0)[1]; double reg2=2.1547005383792514621*e.pos(1)[0]; double reg3=0.15470053837925146212*e.pos(0)[0]; double reg4=2.1547005383792514621*e.pos(1)[1];
   double reg5=2.1547005383792514621*e.pos(0)[0]; double reg6=0.15470053837925146212*e.pos(1)[0]; double reg7=0.15470053837925146212*e.pos(0)[1]; reg1=reg0+reg1; reg4=reg4+reg7;
   double reg8=2.3094010767585029242*e.pos(2)[1]; double reg9=2.3094010767585029242*e.pos(2)[0]; reg5=reg6+reg5; reg2=reg2+reg3; reg4=reg4-reg8;
   reg2=reg2-reg9; double reg10=reg8-reg1; double reg11=reg9-reg5; double reg12=pow(reg4,2); double reg13=pow(reg10,2);
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
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.56758792732771999991*e.pos(0)[1]; double reg3=0.78379396366385999995*e.pos(0)[1]; double reg4=0.56758792732771999991*e.pos(1)[0];
   double reg5=0.56758792732771999991*e.pos(1)[1]; double reg6=0.78379396366385999995*e.pos(0)[0]; double reg7=0.78379396366385999995*e.pos(1)[1]; double reg8=reg7+reg2; double reg9=1.3513818909915799999*e.pos(3)[1];
   double reg10=1.3513818909915799999*e.pos(3)[0]; double reg11=reg0+reg1; double reg12=reg5+reg3; double reg13=0.78379396366385999995*e.pos(1)[2]; double reg14=0.56758792732771999991*e.pos(1)[2];
   double reg15=0.78379396366385999995*e.pos(0)[2]; double reg16=0.56758792732771999991*e.pos(0)[2]; double reg17=reg4+reg6; double reg18=0.63369514596091600003*e.pos(1)[1]; double reg19=2.2673902919218320001*e.pos(0)[1];
   reg8=reg8-reg9; double reg20=reg10-reg17; double reg21=1.78379396366386*e.pos(4)[1]; double reg22=1.3513818909915799999*e.pos(3)[2]; double reg23=reg13+reg16;
   double reg24=reg9-reg12; double reg25=1.78379396366386*e.pos(4)[0]; reg11=reg11-reg10; double reg26=0.63369514596091600003*e.pos(1)[0]; double reg27=2.2673902919218320001*e.pos(0)[0];
   double reg28=reg14+reg15; double reg29=2.2673902919218320001*e.pos(0)[2]; reg20=reg25+reg20; double reg30=reg18+reg19; double reg31=0.63369514596091600003*e.pos(1)[2];
   reg7=reg7-reg3; double reg32=0.43241207267228000009*e.pos(4)[1]; double reg33=0.43241207267228000009*e.pos(4)[0]; double reg34=2.9010854378827480001*e.pos(3)[1]; double reg35=0.63369514596091600003*e.pos(0)[1];
   double reg36=2.2673902919218320001*e.pos(1)[1]; double reg37=0.63369514596091600003*e.pos(0)[0]; reg0=reg0-reg6; double reg38=2.2673902919218320001*e.pos(1)[0]; double reg39=reg22-reg28;
   double reg40=reg26+reg27; reg11=reg11+reg25; double reg41=1.78379396366386*e.pos(5)[0]; double reg42=2.9010854378827480001*e.pos(3)[0]; reg24=reg21+reg24;
   reg8=reg8+reg21; double reg43=1.78379396366386*e.pos(5)[1]; reg23=reg23-reg22; double reg44=1.78379396366386*e.pos(4)[2]; double reg45=reg34-reg30;
   double reg46=reg42-reg40; double reg47=0.366304854039084*e.pos(4)[1]; reg7=reg32+reg7; double reg48=0.366304854039084*e.pos(4)[0]; reg13=reg13-reg15;
   double reg49=0.43241207267228000009*e.pos(5)[1]; reg39=reg44+reg39; reg24=reg24-reg43; reg11=reg11-reg41; reg43=reg8-reg43;
   reg23=reg23+reg44; reg8=1.78379396366386*e.pos(5)[2]; double reg50=reg31+reg29; double reg51=2.9010854378827480001*e.pos(3)[2]; double reg52=0.43241207267228000009*e.pos(5)[0];
   reg0=reg0+reg33; double reg53=0.43241207267228000009*e.pos(4)[2]; reg38=reg38+reg37; reg36=reg36+reg35; double reg54=2.2673902919218320001*e.pos(1)[2];
   double reg55=0.63369514596091600003*e.pos(0)[2]; reg41=reg20-reg41; reg13=reg53+reg13; reg20=pow(reg24,2); reg49=reg7-reg49;
   reg52=reg0-reg52; reg0=pow(reg41,2); reg39=reg39-reg8; reg7=0.43241207267228000009*e.pos(5)[2]; reg8=reg23-reg8;
   reg23=reg37-reg26; reg54=reg54+reg55; double reg56=pow(reg43,2); reg36=reg36-reg34; double reg57=pow(reg11,2);
   double reg58=0.78379396366385999995*e.pos(2)[0]; reg38=reg38-reg42; double reg59=2.710505431213761085e-20*e.pos(3)[0]; double reg60=0.366304854039084*e.pos(4)[2]; double reg61=0.78379396366385999995*e.pos(2)[1];
   double reg62=reg51-reg50; double reg63=reg35-reg18; double reg64=0.366304854039084*e.pos(5)[1]; reg45=reg45+reg47; double reg65=2.710505431213761085e-20*e.pos(3)[1];
   double reg66=0.366304854039084*e.pos(5)[0]; reg46=reg46+reg48; double reg67=pow(reg39,2); reg23=reg23-reg59; reg63=reg63-reg65;
   double reg68=reg55-reg31; double reg69=2.710505431213761085e-20*e.pos(3)[2]; double reg70=3.2673902919218320001*e.pos(4)[1]; double reg71=3.2673902919218320001*e.pos(4)[0]; reg54=reg54-reg51;
   reg36=reg47+reg36; reg38=reg48+reg38; double reg72=0.366304854039084*e.pos(5)[2]; reg62=reg62+reg60; reg45=reg45-reg64;
   double reg73=pow(reg52,2); reg46=reg46-reg66; double reg74=pow(reg49,2); double reg75=0.56758792732771999991*e.pos(2)[1]; reg7=reg13-reg7;
   reg13=0.56758792732771999991*e.pos(2)[0]; double reg76=1.78379396366386*e.pos(3)[1]; double reg77=reg2+reg61; reg56=reg57+reg56; reg57=pow(reg8,2);
   double reg78=reg1+reg58; double reg79=1.78379396366386*e.pos(3)[0]; double reg80=0.78379396366385999995*e.pos(2)[2]; reg20=reg0+reg20; reg0=0.56758792732771999991*e.pos(2)[2];
   double reg81=pow(reg45,2); reg62=reg62-reg72; reg57=reg56+reg57; reg68=reg68-reg69; reg66=reg38-reg66;
   reg38=3.2673902919218320001*e.pos(5)[1]; reg63=reg70+reg63; reg56=0.63369514596091600003*e.pos(2)[0]; double reg82=reg3+reg75; double reg83=3.2673902919218320001*e.pos(5)[0];
   reg23=reg71+reg23; double reg84=3.2673902919218320001*e.pos(4)[2]; reg67=reg20+reg67; reg77=reg77-reg76; reg20=pow(reg7,2);
   reg54=reg60+reg54; reg58=reg58-reg6; double reg85=0.43241207267228000009*e.pos(3)[0]; reg64=reg36-reg64; reg36=1.78379396366386*e.pos(3)[2];
   reg78=reg78-reg79; double reg86=reg16+reg80; reg3=reg61-reg3; reg61=0.43241207267228000009*e.pos(3)[1]; reg74=reg73+reg74;
   reg73=pow(reg46,2); double reg87=0.63369514596091600003*e.pos(2)[1]; reg6=reg6+reg13; reg78=reg25+reg78; double reg88=1.3513818909915799999*e.pos(5)[0];
   double reg89=pow(reg64,2); reg79=reg79+reg6; double reg90=0.366304854039084*e.pos(3)[0]; reg20=reg74+reg20; reg74=pow(reg66,2);
   reg27=reg27+reg56; reg72=reg54-reg72; reg54=pow(reg62,2); reg57=pow(reg57,0.5); reg67=pow(reg67,0.5);
   double reg91=0.63369514596091600003*e.pos(2)[2]; double reg92=0.366304854039084*e.pos(3)[1]; reg86=reg86-reg36; reg19=reg19+reg87; reg85=reg58-reg85;
   reg61=reg3-reg61; reg3=reg15+reg0; reg15=reg80-reg15; reg58=0.43241207267228000009*e.pos(3)[2]; reg81=reg73+reg81;
   reg73=3.2673902919218320001*e.pos(5)[2]; reg68=reg84+reg68; reg80=1.3513818909915799999*e.pos(5)[1]; reg38=reg63-reg38; reg76=reg76+reg82;
   reg83=reg23-reg83; reg77=reg21+reg77; reg25=reg25-reg79; reg36=reg36+reg3; reg23=reg27+reg90;
   reg21=reg21-reg76; reg63=reg19+reg92; double reg93=0.366304854039084*e.pos(3)[2]; reg29=reg29+reg91; reg89=reg74+reg89;
   reg86=reg44+reg86; reg74=2.2673902919218320001*e.pos(2)[1]; reg33=reg85+reg33; reg32=reg61+reg32; reg61=1.3513818909915799999*e.pos(5)[2];
   reg85=2.2673902919218320001*e.pos(2)[0]; reg58=reg15-reg58; reg77=reg77-reg80; reg73=reg68-reg73; reg15=reg24/reg67;
   reg68=pow(reg38,2); double reg94=pow(reg83,2); reg78=reg78-reg88; double reg95=3.2673902919218320001*e.pos(3)[1]; double reg96=reg35-reg87;
   double reg97=3.2673902919218320001*e.pos(3)[0]; double reg98=reg37-reg56; reg54=reg81+reg54; reg81=reg41/reg67; double reg99=pow(reg72,2);
   reg20=pow(reg20,0.5); double reg100=reg11/reg57; double reg101=reg43/reg57; reg86=reg86-reg61; double reg102=reg47-reg63;
   double reg103=reg101*reg77; double reg104=reg48-reg23; double reg105=reg100*reg78; reg54=pow(reg54,0.5); double reg106=2.9010854378827480001*e.pos(5)[0];
   double reg107=2.2673902919218320001*e.pos(2)[2]; reg67=reg39/reg67; reg74=reg35+reg74; double reg108=reg81*reg33; reg85=reg37+reg85;
   double reg109=reg15*reg32; reg53=reg58+reg53; reg58=pow(reg73,2); reg68=reg94+reg68; reg94=3.2673902919218320001*e.pos(3)[2];
   double reg110=reg55-reg91; reg95=reg96-reg95; reg97=reg98-reg97; reg96=2.9010854378827480001*e.pos(5)[1]; reg99=reg89+reg99;
   reg89=reg29+reg93; reg98=reg52/reg20; double reg111=reg49/reg20; reg57=reg8/reg57; reg25=reg88+reg25;
   reg21=reg80+reg21; reg44=reg44-reg36; reg107=reg55+reg107; reg92=reg74-reg92; reg102=reg102+reg96;
   reg90=reg85-reg90; reg74=reg60-reg89; reg85=2.9010854378827480001*e.pos(5)[2]; reg103=reg105+reg103; reg58=reg68+reg58;
   reg104=reg106+reg104; reg94=reg110-reg94; reg68=5.42101086242752217e-20*e.pos(5)[1]; reg70=reg95+reg70; reg71=reg97+reg71;
   reg95=5.42101086242752217e-20*e.pos(5)[0]; reg97=reg45/reg54; reg99=pow(reg99,0.5); reg105=reg46/reg54; reg110=reg57*reg86;
   reg109=reg108+reg109; reg20=reg7/reg20; reg108=reg98*reg25; double reg112=reg111*reg21; reg44=reg61+reg44;
   double reg113=reg67*reg53; reg93=reg107-reg93; reg54=reg62/reg54; reg92=reg47+reg92; reg90=reg48+reg90;
   reg47=reg105*reg104; reg48=reg97*reg102; reg107=reg20*reg44; reg74=reg74+reg85; reg58=pow(reg58,0.5);
   reg113=reg109+reg113; reg112=reg108+reg112; reg108=reg66/reg99; reg109=reg64/reg99; double reg114=5.42101086242752217e-20*e.pos(5)[2];
   reg84=reg94+reg84; reg70=reg70-reg68; reg71=reg71-reg95; reg110=reg103+reg110; reg94=reg108*reg71;
   reg103=reg100*reg110; double reg115=reg109*reg70; reg99=reg72/reg99; reg84=reg84-reg114; double reg116=reg101*reg110;
   double reg117=reg15*reg113; double reg118=reg81*reg113; double reg119=reg83/reg58; double reg120=reg38/reg58; double reg121=reg54*reg74;
   reg107=reg112+reg107; reg48=reg47+reg48; reg93=reg60+reg93; reg92=reg92-reg96; reg90=reg90-reg106;
   reg47=reg99*reg84; reg60=reg120*reg92; reg116=reg77-reg116; reg121=reg48+reg121; reg58=reg73/reg58;
   reg48=reg119*reg90; reg118=reg33-reg118; reg112=reg98*reg107; double reg122=reg67*reg113; double reg123=reg111*reg107;
   reg103=reg78-reg103; double reg124=reg57*reg110; reg117=reg32-reg117; reg115=reg94+reg115; reg93=reg93-reg85;
   reg94=pow(reg118,2); reg122=reg53-reg122; reg124=reg86-reg124; double reg125=pow(reg116,2); double reg126=pow(reg117,2);
   reg47=reg115+reg47; reg115=reg97*reg121; double reg127=reg20*reg107; reg123=reg21-reg123; double reg128=reg105*reg121;
   double reg129=reg58*reg93; reg60=reg48+reg60; reg112=reg25-reg112; reg48=pow(reg103,2); reg94=reg126+reg94;
   reg127=reg44-reg127; reg125=reg48+reg125; reg48=pow(reg123,2); reg126=pow(reg124,2); double reg130=pow(reg122,2);
   double reg131=reg109*reg47; double reg132=pow(reg112,2); reg115=reg102-reg115; double reg133=reg54*reg121; reg128=reg104-reg128;
   double reg134=reg108*reg47; reg129=reg60+reg129; reg134=reg71-reg134; reg131=reg70-reg131; reg60=reg99*reg47;
   double reg135=reg120*reg129; double reg136=reg119*reg129; reg94=reg130+reg94; reg130=pow(reg128,2); reg48=reg132+reg48;
   reg133=reg74-reg133; reg126=reg125+reg126; reg125=pow(reg115,2); reg132=pow(reg127,2); double reg137=pow(reg131,2);
   reg60=reg84-reg60; double reg138=pow(reg134,2); reg132=reg48+reg132; reg48=reg58*reg129; reg136=reg90-reg136;
   reg135=reg92-reg135; reg125=reg130+reg125; reg130=pow(reg133,2); reg94=pow(reg94,0.5); reg126=pow(reg126,0.5);
   reg132=pow(reg132,0.5); double reg139=pow(reg136,2); double reg140=pow(reg135,2); reg130=reg125+reg130; reg48=reg93-reg48;
   reg118=reg118/reg94; reg116=reg116/reg126; reg137=reg138+reg137; reg117=reg117/reg94; reg103=reg103/reg126;
   reg125=pow(reg60,2); reg138=pow(reg48,2); reg140=reg139+reg140; reg33=reg33*reg118; reg32=reg32*reg117;
   reg118=reg41*reg118; reg117=reg24*reg117; reg15=reg24*reg15; reg81=reg41*reg81; reg125=reg137+reg125;
   reg94=reg122/reg94; reg126=reg124/reg126; reg78=reg78*reg103; reg77=reg77*reg116; reg100=reg11*reg100;
   reg103=reg11*reg103; reg116=reg43*reg116; reg101=reg43*reg101; reg130=pow(reg130,0.5); reg112=reg112/reg132;
   reg123=reg123/reg132; reg128=reg128/reg130; reg125=pow(reg125,0.5); reg15=reg81+reg15; reg115=reg115/reg130;
   reg111=reg49*reg111; reg67=reg39*reg67; reg98=reg52*reg98; reg25=reg25*reg112; reg77=reg78+reg77;
   reg21=reg21*reg123; reg86=reg86*reg126; reg53=reg53*reg94; reg57=reg8*reg57; reg123=reg49*reg123;
   reg112=reg52*reg112; reg101=reg100+reg101; reg138=reg140+reg138; reg126=reg8*reg126; reg116=reg103+reg116;
   reg94=reg39*reg94; reg117=reg118+reg117; reg32=reg33+reg32; reg132=reg127/reg132; reg104=reg104*reg128;
   reg57=reg101+reg57; reg138=pow(reg138,0.5); reg130=reg133/reg130; reg128=reg46*reg128; reg53=reg32+reg53;
   reg102=reg102*reg115; reg126=reg116+reg126; reg8=reg7*reg132; reg123=reg112+reg123; reg94=reg117+reg94;
   reg105=reg46*reg105; reg132=reg44*reg132; reg97=reg45*reg97; reg115=reg45*reg115; reg86=reg77+reg86;
   reg111=reg98+reg111; reg20=reg7*reg20; reg134=reg134/reg125; reg131=reg131/reg125; reg67=reg15+reg67;
   reg21=reg25+reg21; reg115=reg128+reg115; reg7=reg62*reg130; reg130=reg74*reg130; reg102=reg104+reg102;
   reg54=reg62*reg54; reg97=reg105+reg97; reg132=reg21+reg132; reg8=reg123+reg8; reg125=reg60/reg125;
   reg126=reg110*reg126; reg70=reg70*reg131; reg11=reg66*reg134; reg20=reg111+reg20; reg131=reg64*reg131;
   reg94=reg113*reg94; reg53=reg67*reg53; reg136=reg136/reg138; reg86=reg57*reg86; reg134=reg71*reg134;
   reg135=reg135/reg138; reg108=reg66*reg108; reg109=reg64*reg109; reg84=reg84*reg125; reg70=reg134+reg70;
   reg54=reg97+reg54; reg109=reg108+reg109; reg99=reg72*reg99; reg132=reg20*reg132; reg94=reg53-reg94;
   reg8=reg107*reg8; reg7=reg115+reg7; reg15=reg38*reg135; reg20=reg83*reg136; reg130=reg102+reg130;
   reg126=reg86-reg126; reg136=reg90*reg136; reg119=reg83*reg119; reg135=reg92*reg135; reg120=reg38*reg120;
   reg125=reg72*reg125; reg131=reg11+reg131; reg138=reg48/reg138; reg93=reg93*reg138; reg135=reg136+reg135;
   reg11=0.009463616120767210603*reg94; reg21=0.088847818743090689935*reg94; reg24=0.088847818743090689935*reg126; reg15=reg20+reg15; reg20=0.005384432036113586778*reg126;
   reg25=0.021537728144454347112*reg94; reg32=0.021537728144454347112*reg126; reg138=reg73*reg138; reg58=reg73*reg58; reg120=reg119+reg120;
   reg33=0.009463616120767210603*reg126; reg7=reg121*reg7; reg130=reg54*reg130; reg8=reg132-reg8; reg125=reg131+reg125;
   reg38=0.005384432036113586778*reg94; reg84=reg70+reg84; reg99=reg109+reg99; reg21=reg32+reg21; reg11=reg20+reg11;
   reg20=reg38+reg20; reg39=0.009463616120767210603*reg8; reg41=0.021537728144454347112*reg8; reg24=reg25+reg24; reg43=0.005384432036113586778*reg8;
   reg25=reg32+reg25; reg38=reg33+reg38; reg32=0.088847818743090689935*reg8; reg58=reg120+reg58; reg138=reg15+reg138;
   reg7=reg130-reg7; reg125=reg47*reg125; reg84=reg99*reg84; reg93=reg135+reg93; reg15=0.016449618187943419918*reg7;
   reg32=reg25+reg32; reg125=reg84-reg125; reg38=reg38+reg43; reg93=reg58*reg93; reg39=reg20+reg39;
   reg24=reg24+reg41; reg20=0.028457289286966203713*reg7; reg25=0.0018441552587796664112*reg7; reg33=0.0041124045469858549794*reg7; reg11=reg43+reg11;
   reg138=reg129*reg138; reg21=reg41+reg21; reg21=reg15+reg21; reg25=reg24+reg25; reg24=0.0018441552587796664109*reg125;
   reg41=0.016449618187943419918*reg125; reg43=0.016449618187943419916*reg125; reg44=0.004112404546985854979*reg125; reg38=reg20-reg38; reg138=reg93-reg138;
   reg20=0.028457289286966203713*reg125; reg11=reg11+reg33; reg39=reg33+reg39; reg33=0.0041124045469858549794*reg125; reg15=reg32+reg15;
   reg32=0.028457289286966203713*reg138; reg33=reg39+reg33; reg11=reg20-reg11; reg24=reg21+reg24; reg20=0.0018441552587796664111*reg138;
   reg41=reg25+reg41; reg21=0.016449618187943419918*reg138; reg44=reg38-reg44; reg25=0.0041124045469858549794*reg138; reg43=reg15+reg43;
   Ne(0,15)+=reg21+reg24; Ne(1,16)+=reg21+reg24; Ne(2,17)+=reg21+reg24; Ne(0,6)+=reg32-reg33; Ne(1,7)+=reg32-reg33;
   Ne(2,8)+=reg32-reg33; Ne(0,0)+=reg44-reg25; Ne(1,1)+=reg44-reg25; Ne(2,2)+=reg44-reg25; Ne(0,12)+=reg41+reg21;
   Ne(1,13)+=reg41+reg21; Ne(2,14)+=reg41+reg21; Ne(0,9)+=reg43+reg20; Ne(1,10)+=reg43+reg20; Ne(2,11)+=reg43+reg20;
   Ne(0,3)+=reg11-reg25; Ne(1,4)+=reg11-reg25; Ne(2,5)+=reg11-reg25;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(1)[0]; double reg1=0.12200846792814621817*e.pos(1)[0]; double reg2=0.12200846792814621817*e.pos(1)[1]; double reg3=0.36602540378443865451*e.pos(0)[0]; double reg4=0.36602540378443865451*e.pos(0)[1];
   double reg5=0.12200846792814621817*e.pos(0)[0]; double reg6=0.12200846792814621817*e.pos(0)[1]; double reg7=0.36602540378443865451*e.pos(1)[1]; double reg8=1.3660254037844385386*e.pos(2)[0]; double reg9=reg1+reg3;
   double reg10=0.45534180126147951289*e.pos(2)[1]; double reg11=reg7+reg6; double reg12=reg0+reg5; double reg13=0.45534180126147951289*e.pos(2)[0]; double reg14=0.12200846792814621817*e.pos(1)[2];
   double reg15=0.36602540378443865451*e.pos(0)[2]; double reg16=1.3660254037844385386*e.pos(2)[1]; double reg17=1.3660254037844385386*e.pos(1)[1]; double reg18=reg2+reg4; double reg19=1.3660254037844385386*e.pos(1)[0];
   double reg20=0.45534180126147951289*e.pos(0)[1]; double reg21=0.45534180126147951289*e.pos(0)[0]; double reg22=0.12200846792814621817*e.pos(0)[2]; double reg23=0.36602540378443865451*e.pos(1)[2]; double reg24=1.3660254037844385386*e.pos(2)[2];
   double reg25=0.45534180126147951289*e.pos(1)[0]; double reg26=0.45534180126147951289*e.pos(1)[1]; double reg27=reg14+reg15; double reg28=reg12+reg13; double reg29=1.3660254037844385386*e.pos(3)[0];
   double reg30=reg11+reg10; double reg31=1.3660254037844385386*e.pos(3)[1]; double reg32=reg23+reg22; double reg33=0.45534180126147951289*e.pos(2)[2]; double reg34=0.12200846792814621817*e.pos(2)[0];
   double reg35=0.12200846792814621817*e.pos(2)[1]; double reg36=0.45534180126147951289*e.pos(0)[2]; double reg37=reg21+reg19; double reg38=reg20+reg17; double reg39=1.3660254037844385386*e.pos(1)[2];
   double reg40=1.3660254037844385386*e.pos(0)[1]; double reg41=1.3660254037844385386*e.pos(0)[0]; double reg42=0.45534180126147951289*e.pos(3)[0]; reg9=reg9+reg8; reg18=reg18+reg16;
   double reg43=0.45534180126147951289*e.pos(3)[1]; double reg44=1.3660254037844385386*e.pos(0)[2]; double reg45=0.45534180126147951289*e.pos(1)[2]; double reg46=0.36602540378443865451*e.pos(2)[1]; reg37=reg34+reg37;
   double reg47=0.36602540378443865451*e.pos(3)[0]; reg38=reg35+reg38; double reg48=0.36602540378443865451*e.pos(3)[1]; double reg49=reg36+reg39; double reg50=0.12200846792814621817*e.pos(2)[2];
   double reg51=1.3660254037844385386*e.pos(3)[2]; double reg52=reg32+reg33; double reg53=reg30+reg31; reg9=reg9+reg42; double reg54=0.48803387171258487271*e.pos(4)[0];
   double reg55=reg28+reg29; double reg56=0.36602540378443865451*e.pos(2)[0]; double reg57=reg25+reg41; reg18=reg18+reg43; reg27=reg27+reg24;
   double reg58=0.45534180126147951289*e.pos(3)[2]; double reg59=reg26+reg40; double reg60=0.48803387171258487271*e.pos(4)[1]; double reg61=reg59+reg46; reg27=reg27+reg58;
   double reg62=0.48803387171258487271*e.pos(4)[2]; double reg63=reg57+reg56; reg37=reg37+reg47; double reg64=1.8213672050459180516*e.pos(4)[0]; double reg65=reg52+reg51;
   double reg66=0.12200846792814621817*e.pos(3)[1]; double reg67=0.36602540378443865451*e.pos(2)[2]; double reg68=reg45+reg44; reg9=reg9-reg54; double reg69=0.66666666666666670528*e.pos(5)[0];
   double reg70=reg60-reg53; reg38=reg38+reg48; double reg71=1.8213672050459180516*e.pos(4)[1]; reg18=reg18-reg60; double reg72=0.12200846792814621817*e.pos(3)[0];
   double reg73=0.66666666666666670528*e.pos(5)[1]; double reg74=0.36602540378443865451*e.pos(3)[2]; double reg75=reg54-reg55; reg49=reg50+reg49; double reg76=reg68+reg67;
   reg9=reg9+reg69; double reg77=1.8213672050459180516*e.pos(6)[0]; reg75=reg69+reg75; reg18=reg18+reg73; double reg78=reg62-reg65;
   double reg79=1.8213672050459180516*e.pos(6)[1]; reg70=reg73+reg70; reg38=reg38-reg71; reg37=reg37-reg64; double reg80=reg72+reg63;
   reg49=reg49+reg74; double reg81=reg66+reg61; double reg82=1.8213672050459180516*e.pos(4)[2]; double reg83=0.66666666666666670528*e.pos(5)[2]; reg27=reg27-reg62;
   double reg84=0.12200846792814621817*e.pos(3)[2]; double reg85=0.48803387171258487271*e.pos(6)[1]; reg38=reg73+reg38; reg18=reg18-reg79; double reg86=0.66666666666666670528*e.pos(7)[1];
   reg75=reg77+reg75; double reg87=reg3+reg25; double reg88=reg64-reg80; reg70=reg79+reg70; reg9=reg9-reg77;
   double reg89=0.66666666666666670528*e.pos(7)[0]; double reg90=reg4+reg26; double reg91=0.48803387171258487271*e.pos(6)[0]; reg37=reg69+reg37; reg78=reg83+reg78;
   double reg92=1.8213672050459180516*e.pos(6)[2]; double reg93=reg0+reg21; double reg94=reg84+reg76; reg27=reg27+reg83; double reg95=reg7+reg20;
   reg49=reg49-reg82; double reg96=reg71-reg81; reg19=reg5+reg19; reg38=reg38-reg85; reg5=0.48803387171258487271*e.pos(6)[2];
   reg75=reg75-reg89; reg49=reg83+reg49; reg17=reg6+reg17; reg88=reg69+reg88; reg70=reg70-reg86;
   reg37=reg37-reg91; reg78=reg92+reg78; reg34=reg34+reg93; reg35=reg35+reg95; reg6=reg23+reg36;
   reg69=reg82-reg94; reg96=reg73+reg96; reg73=reg15+reg45; reg16=reg16+reg90; double reg97=0.66666666666666670528*e.pos(7)[2];
   reg9=reg9-reg89; reg27=reg27-reg92; reg8=reg8+reg87; reg18=reg18-reg86; reg72=reg8+reg72;
   reg38=reg38-reg86; reg8=0.66666666666666670528*e.pos(4)[0]; reg41=reg1+reg41; reg37=reg37-reg89; reg49=reg49-reg5;
   reg34=reg29+reg34; reg66=reg16+reg66; reg1=0.66666666666666670528*e.pos(4)[1]; reg78=reg78-reg97; reg35=reg31+reg35;
   reg88=reg91+reg88; reg19=reg13+reg19; reg16=pow(reg18,2); reg29=pow(reg9,2); reg31=pow(reg75,2);
   reg40=reg2+reg40; reg69=reg83+reg69; reg27=reg27-reg97; reg50=reg50+reg6; reg39=reg22+reg39;
   reg24=reg24+reg73; reg17=reg10+reg17; reg96=reg85+reg96; reg2=pow(reg70,2); reg22=1.8213672050459180516*e.pos(5)[1];
   reg50=reg51+reg50; reg51=pow(reg37,2); reg83=1.8213672050459180516*e.pos(5)[0]; reg66=reg66-reg1; double reg98=0.48803387171258487271*e.pos(5)[0];
   double reg99=0.66666666666666670528*e.pos(4)[2]; reg35=reg35-reg1; reg84=reg24+reg84; reg24=0.48803387171258487271*e.pos(5)[1]; reg69=reg5+reg69;
   reg34=reg34-reg8; reg86=reg96-reg86; reg44=reg14+reg44; reg2=reg31+reg2; reg14=pow(reg27,2);
   reg31=reg48+reg17; reg96=reg47+reg19; reg39=reg33+reg39; reg89=reg88-reg89; reg88=pow(reg78,2);
   reg49=reg49-reg97; reg40=reg46+reg40; reg16=reg29+reg16; reg41=reg56+reg41; reg29=pow(reg38,2);
   reg72=reg72-reg8; double reg100=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg99; double reg101=reg74+reg39; double reg102=reg1+reg31;
   reg97=reg69-reg97; reg69=pow(reg86,2); double reg103=reg8+reg96; double reg104=0.66666666666666670528*e.pos(6)[1]; reg66=reg66-reg22;
   reg14=reg16+reg14; reg16=pow(reg89,2); double reg105=pow(reg49,2); double reg106=0.66666666666666670528*e.pos(6)[0]; reg29=reg51+reg29;
   reg72=reg72-reg83; reg50=reg50-reg99; reg51=0.48803387171258487271*e.pos(5)[2]; reg44=reg67+reg44; reg35=reg35-reg24;
   reg34=reg34-reg98; double reg107=reg42+reg41; reg88=reg2+reg88; reg2=reg43+reg40; reg88=pow(reg88,0.5);
   reg69=reg16+reg69; reg16=1.8213672050459180516*e.pos(7)[1]; double reg108=reg83-reg103; double reg109=reg58+reg44; double reg110=0.48803387171258487271*e.pos(7)[1];
   reg66=reg66+reg104; reg35=reg104+reg35; reg14=pow(reg14,0.5); double reg111=1.8213672050459180516*e.pos(7)[0]; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg72=reg106+reg72; reg34=reg106+reg34; reg105=reg29+reg105; reg29=reg22-reg102; double reg113=reg99+reg101;
   double reg114=pow(reg97,2); reg1=reg1+reg2; reg50=reg50-reg51; double reg115=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg100;
   reg8=reg8+reg107; double reg116=reg9/reg14; double reg117=reg75/reg88; double reg118=reg18/reg14; double reg119=reg98-reg8;
   double reg120=reg100-reg113; double reg121=reg70/reg88; double reg122=reg24-reg1; reg72=reg72-reg112; reg34=reg34-reg111;
   reg105=pow(reg105,0.5); reg84=reg84+reg115; double reg123=0.48803387171258487271*e.pos(7)[2]; reg108=reg106+reg108; reg99=reg99+reg109;
   reg29=reg104+reg29; reg114=reg69+reg114; reg66=reg66-reg110; reg35=reg35-reg16; reg69=1.8213672050459180516*e.pos(7)[2];
   reg50=reg115+reg50; reg14=reg27/reg14; reg108=reg112+reg108; reg119=reg106+reg119; reg120=reg115+reg120;
   reg122=reg104+reg122; reg104=reg51-reg99; reg106=reg38/reg105; double reg124=reg37/reg105; reg29=reg110+reg29;
   reg50=reg50-reg69; reg114=pow(reg114,0.5); reg84=reg84-reg123; double reg125=reg121*reg35; double reg126=reg118*reg66;
   double reg127=reg117*reg34; double reg128=reg116*reg72; reg88=reg78/reg88; double reg129=reg106*reg29; reg120=reg123+reg120;
   double reg130=reg124*reg108; reg125=reg127+reg125; reg126=reg128+reg126; reg122=reg16+reg122; reg127=reg14*reg84;
   reg104=reg115+reg104; reg115=reg89/reg114; reg105=reg49/reg105; reg128=reg86/reg114; double reg131=reg88*reg50;
   reg119=reg111+reg119; double reg132=reg105*reg120; double reg133=reg128*reg122; reg131=reg125+reg131; reg129=reg130+reg129;
   reg114=reg97/reg114; reg125=reg115*reg119; reg127=reg126+reg127; reg104=reg69+reg104; reg126=reg118*reg127;
   reg130=reg116*reg127; reg132=reg129+reg132; reg133=reg125+reg133; reg125=reg114*reg104; reg129=reg121*reg131;
   double reg134=reg117*reg131; reg129=reg35-reg129; double reg135=reg124*reg132; reg134=reg34-reg134; double reg136=reg106*reg132;
   reg125=reg133+reg125; reg133=reg88*reg131; reg126=reg66-reg126; reg130=reg72-reg130; double reg137=reg14*reg127;
   double reg138=pow(reg126,2); double reg139=pow(reg130,2); double reg140=reg105*reg132; reg136=reg29-reg136; reg137=reg84-reg137;
   double reg141=reg128*reg125; reg135=reg108-reg135; double reg142=reg115*reg125; reg133=reg50-reg133; double reg143=pow(reg134,2);
   double reg144=pow(reg129,2); reg142=reg119-reg142; double reg145=reg114*reg125; reg141=reg122-reg141; double reg146=pow(reg135,2);
   double reg147=pow(reg136,2); reg140=reg120-reg140; double reg148=pow(reg137,2); reg138=reg139+reg138; reg139=pow(reg133,2);
   reg144=reg143+reg144; reg147=reg146+reg147; reg143=pow(reg140,2); reg139=reg144+reg139; reg148=reg138+reg148;
   reg145=reg104-reg145; reg138=pow(reg142,2); reg144=pow(reg141,2); reg148=pow(reg148,0.5); reg144=reg138+reg144;
   reg139=pow(reg139,0.5); reg143=reg147+reg143; reg138=pow(reg145,2); reg130=reg130/reg148; reg126=reg126/reg148;
   reg134=reg134/reg139; reg143=pow(reg143,0.5); reg129=reg129/reg139; reg138=reg144+reg138; reg34=reg34*reg134;
   reg35=reg35*reg129; reg139=reg133/reg139; reg121=reg70*reg121; reg117=reg75*reg117; reg133=reg18*reg126;
   reg144=reg9*reg130; reg148=reg137/reg148; reg126=reg66*reg126; reg130=reg72*reg130; reg118=reg18*reg118;
   reg116=reg9*reg116; reg138=pow(reg138,0.5); reg135=reg135/reg143; reg136=reg136/reg143; reg129=reg70*reg129;
   reg134=reg75*reg134; reg84=reg84*reg148; reg126=reg130+reg126; reg133=reg144+reg133; reg148=reg27*reg148;
   reg9=reg38*reg136; reg18=reg37*reg135; reg14=reg27*reg14; reg118=reg116+reg118; reg106=reg38*reg106;
   reg121=reg117+reg121; reg143=reg140/reg143; reg88=reg78*reg88; reg136=reg29*reg136; reg50=reg50*reg139;
   reg35=reg34+reg35; reg129=reg134+reg129; reg139=reg78*reg139; reg141=reg141/reg138; reg124=reg37*reg124;
   reg142=reg142/reg138; reg135=reg108*reg135; reg119=reg119*reg142; reg120=reg120*reg143; reg122=reg122*reg141;
   reg138=reg145/reg138; reg14=reg118+reg14; reg142=reg89*reg142; reg136=reg135+reg136; reg9=reg18+reg9;
   reg143=reg49*reg143; reg141=reg86*reg141; reg50=reg35+reg50; reg139=reg129+reg139; reg115=reg89*reg115;
   reg88=reg121+reg88; reg128=reg86*reg128; reg106=reg124+reg106; reg148=reg133+reg148; reg105=reg49*reg105;
   reg84=reg126+reg84; reg122=reg119+reg122; reg139=reg131*reg139; reg84=reg14*reg84; reg104=reg104*reg138;
   reg50=reg88*reg50; reg128=reg115+reg128; reg114=reg97*reg114; reg141=reg142+reg141; reg138=reg97*reg138;
   reg120=reg136+reg120; reg148=reg127*reg148; reg143=reg9+reg143; reg105=reg106+reg105; reg138=reg141+reg138;
   reg104=reg122+reg104; reg139=reg50-reg139; reg148=reg84-reg148; reg114=reg128+reg114; reg143=reg132*reg143;
   reg120=reg105*reg120; reg9=0.024056261216234395431*reg148; reg14=0.13144585576580215187*reg139; reg18=0.024056261216234395431*reg139; reg27=0.024056261216234409915*reg139;
   reg29=0.13144585576580215187*reg148; reg34=0.04166666666666666908*reg148; reg35=0.035220810900864524453*reg148; reg37=0.024056261216234409915*reg148; reg38=0.035220810900864524453*reg139;
   reg138=reg125*reg138; reg104=reg114*reg104; reg143=reg120-reg143; reg49=0.04166666666666666908*reg139; reg50=reg35+reg38;
   reg66=0.024056261216234409915*reg143; reg18=reg18-reg34; reg70=0.13144585576580215187*reg143; reg38=reg38+reg29; reg9=reg9-reg49;
   reg29=reg29+reg14; reg72=0.035220810900864524453*reg143; reg27=reg34+reg27; reg34=0.024056261216234395431*reg143; reg14=reg35+reg14;
   reg138=reg104-reg138; reg35=0.04166666666666666908*reg143; reg49=reg37+reg49; reg37=0.035220810900864524453*reg138; reg38=reg70+reg38;
   reg75=0.13144585576580215187*reg138; reg29=reg29+reg72; reg14=reg72+reg14; reg70=reg50+reg70; reg66=reg18-reg66;
   reg18=0.024056261216234409915*reg138; reg9=reg9-reg35; reg50=0.04166666666666666908*reg138; reg27=reg34-reg27; reg34=0.024056261216234395431*reg138;
   reg35=reg49+reg35; Ne(0,9)+=reg66-reg50; Ne(1,10)+=reg66-reg50; Ne(2,11)+=reg66-reg50; Ne(0,12)+=reg70+reg75;
   Ne(1,13)+=reg70+reg75; Ne(2,14)+=reg70+reg75; Ne(0,6)+=reg9-reg18; Ne(1,7)+=reg9-reg18; Ne(2,8)+=reg9-reg18;
   Ne(0,15)+=reg38+reg37; Ne(1,16)+=reg38+reg37; Ne(2,17)+=reg38+reg37; Ne(0,3)+=reg27-reg50; Ne(1,4)+=reg27-reg50;
   Ne(2,5)+=reg27-reg50; Ne(0,0)+=reg34-reg35; Ne(1,1)+=reg34-reg35; Ne(2,2)+=reg34-reg35; Ne(0,18)+=reg37+reg29;
   Ne(1,19)+=reg37+reg29; Ne(2,20)+=reg37+reg29; Ne(0,21)+=reg75+reg14; Ne(1,22)+=reg75+reg14; Ne(2,23)+=reg75+reg14;

}

};
} // namespace LMT

