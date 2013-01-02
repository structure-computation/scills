
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
   double reg5=pow(reg4,2); reg2=reg3+reg2; reg5=reg2+reg5; reg5=pow(reg5,0.5); reg2=reg0/reg5;
   reg3=reg1/reg5; double reg6=e.pos(2)[0]-e.pos(0)[0]; double reg7=e.pos(2)[1]-e.pos(0)[1]; double reg8=reg2*reg6; double reg9=reg3*reg7;
   double reg10=e.pos(2)[2]-e.pos(0)[2]; reg5=reg4/reg5; reg9=reg8+reg9; reg8=reg5*reg10; reg8=reg9+reg8;
   reg9=reg2*reg8; double reg11=reg3*reg8; reg9=reg6-reg9; reg11=reg7-reg11; double reg12=reg5*reg8;
   reg12=reg10-reg12; double reg13=pow(reg9,2); double reg14=pow(reg11,2); reg14=reg13+reg14; reg13=pow(reg12,2);
   reg13=reg14+reg13; reg13=pow(reg13,0.5); reg9=reg9/reg13; reg11=reg11/reg13; reg14=reg0*reg9;
   reg7=reg7*reg11; reg9=reg6*reg9; reg11=reg1*reg11; reg3=reg1*reg3; reg13=reg12/reg13;
   reg2=reg0*reg2; reg5=reg4*reg5; reg10=reg10*reg13; reg7=reg9+reg7; reg13=reg4*reg13;
   reg11=reg14+reg11; reg3=reg2+reg3; reg13=reg11+reg13; reg5=reg3+reg5; reg10=reg7+reg10;
   reg0=reg5*reg10; reg1=reg8*reg13; reg1=reg0-reg1; Ne(0,0)+=0.16666666666666668517*reg1; Ne(1,1)+=0.16666666666666668517*reg1;
   Ne(2,2)+=0.16666666666666668517*reg1; Ne(0,3)+=0.16666666666666665741*reg1; Ne(0,6)+=0.16666666666666665741*reg1; Ne(1,4)+=0.16666666666666665741*reg1; Ne(1,7)+=0.16666666666666665741*reg1;
   Ne(2,5)+=0.16666666666666665741*reg1; Ne(2,8)+=0.16666666666666665741*reg1;

}
 /// pour les Quad
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.21132486540518713447*e.pos(0)[0]; double reg1=0.21132486540518713447*e.pos(1)[0]; double reg2=0.21132486540518713447*e.pos(0)[1]; double reg3=0.21132486540518713447*e.pos(1)[1]; double reg4=0.78867513459481286553*e.pos(1)[1];
   double reg5=0.21132486540518713447*e.pos(1)[2]; double reg6=0.78867513459481286553*e.pos(0)[0]; double reg7=reg3-reg2; double reg8=0.78867513459481286553*e.pos(2)[1]; double reg9=0.21132486540518713447*e.pos(0)[2];
   double reg10=0.78867513459481286553*e.pos(2)[0]; double reg11=reg1-reg0; double reg12=0.78867513459481286553*e.pos(0)[1]; double reg13=0.78867513459481286553*e.pos(1)[0]; double reg14=0.21132486540518713447*e.pos(2)[1];
   double reg15=0.78867513459481286553*e.pos(2)[2]; double reg16=reg13-reg6; double reg17=reg5-reg9; double reg18=0.21132486540518713447*e.pos(2)[0]; double reg19=0.78867513459481286553*e.pos(1)[2];
   double reg20=0.78867513459481286553*e.pos(3)[1]; reg7=reg7+reg8; double reg21=0.78867513459481286553*e.pos(3)[0]; reg11=reg11+reg10; double reg22=reg4-reg12;
   double reg23=0.78867513459481286553*e.pos(0)[2]; double reg24=0.21132486540518713447*e.pos(3)[1]; double reg25=reg19-reg23; double reg26=0.21132486540518713447*e.pos(2)[2]; double reg27=0.21132486540518713447*e.pos(3)[0];
   reg22=reg22+reg14; reg16=reg18+reg16; reg7=reg7-reg20; reg11=reg11-reg21; double reg28=0.78867513459481286553*e.pos(3)[2];
   reg17=reg17+reg15; reg16=reg16-reg27; double reg29=0.21132486540518713447*e.pos(3)[2]; double reg30=pow(reg11,2); reg25=reg26+reg25;
   reg22=reg22-reg24; reg17=reg17-reg28; double reg31=pow(reg7,2); double reg32=pow(reg17,2); double reg33=pow(reg22,2);
   reg31=reg30+reg31; reg25=reg25-reg29; reg30=pow(reg16,2); reg4=reg2+reg4; reg6=reg1+reg6;
   reg33=reg30+reg33; reg13=reg0+reg13; reg32=reg31+reg32; reg12=reg3+reg12; reg0=pow(reg25,2);
   reg0=reg33+reg0; reg19=reg9+reg19; reg6=reg18-reg6; reg23=reg5+reg23; reg4=reg8-reg4;
   reg13=reg10-reg13; reg12=reg14-reg12; reg32=pow(reg32,0.5); reg27=reg13+reg27; reg23=reg26-reg23;
   reg0=pow(reg0,0.5); reg6=reg21+reg6; reg12=reg20+reg12; reg1=reg11/reg32; reg2=reg7/reg32;
   reg24=reg4+reg24; reg19=reg15-reg19; reg23=reg28+reg23; reg32=reg17/reg32; reg3=reg1*reg27;
   reg4=reg1*reg6; reg5=reg2*reg24; reg8=reg2*reg12; reg29=reg19+reg29; reg9=reg16/reg0;
   reg10=reg22/reg0; reg0=reg25/reg0; reg5=reg3+reg5; reg3=reg32*reg29; reg13=reg27*reg9;
   reg14=reg24*reg10; reg8=reg4+reg8; reg4=reg32*reg23; reg15=reg6*reg9; reg14=reg13+reg14;
   reg13=reg29*reg0; reg3=reg5+reg3; reg5=reg12*reg10; reg4=reg8+reg4; reg8=reg2*reg4;
   reg18=reg1*reg4; reg19=reg1*reg3; reg13=reg14+reg13; reg14=reg23*reg0; reg20=reg2*reg3;
   reg5=reg15+reg5; reg15=reg32*reg3; reg21=reg32*reg4; reg14=reg5+reg14; reg8=reg12-reg8;
   reg18=reg6-reg18; reg20=reg24-reg20; reg19=reg27-reg19; reg5=reg10*reg13; reg26=reg9*reg13;
   reg28=pow(reg19,2); reg30=pow(reg20,2); reg31=reg9*reg14; reg33=reg10*reg14; reg15=reg29-reg15;
   reg21=reg23-reg21; double reg34=pow(reg8,2); double reg35=pow(reg18,2); double reg36=reg0*reg13; reg5=reg24-reg5;
   reg26=reg27-reg26; reg36=reg29-reg36; double reg37=pow(reg5,2); double reg38=pow(reg21,2); double reg39=pow(reg26,2);
   reg31=reg6-reg31; reg33=reg12-reg33; double reg40=reg0*reg14; reg34=reg35+reg34; reg35=pow(reg15,2);
   reg30=reg28+reg30; reg40=reg23-reg40; reg28=pow(reg33,2); reg37=reg39+reg37; reg35=reg30+reg35;
   reg38=reg34+reg38; reg30=pow(reg31,2); reg34=pow(reg36,2); reg28=reg30+reg28; reg30=pow(reg40,2);
   reg38=pow(reg38,0.5); reg34=reg37+reg34; reg35=pow(reg35,0.5); reg34=pow(reg34,0.5); reg8=reg8/reg38;
   reg18=reg18/reg38; reg20=reg20/reg35; reg30=reg28+reg30; reg19=reg19/reg35; reg28=reg7*reg20;
   reg2=reg7*reg2; reg37=reg11*reg19; reg19=reg27*reg19; reg35=reg15/reg35; reg20=reg24*reg20;
   reg1=reg11*reg1; reg26=reg26/reg34; reg15=reg6*reg18; reg39=reg12*reg8; reg5=reg5/reg34;
   reg38=reg21/reg38; reg30=pow(reg30,0.5); reg18=reg11*reg18; reg8=reg7*reg8; reg27=reg27*reg26;
   reg24=reg24*reg5; reg34=reg36/reg34; reg26=reg16*reg26; reg5=reg22*reg5; reg2=reg1+reg2;
   reg20=reg19+reg20; reg8=reg18+reg8; reg1=reg17*reg38; reg33=reg33/reg30; reg31=reg31/reg30;
   reg38=reg23*reg38; reg39=reg15+reg39; reg9=reg16*reg9; reg10=reg22*reg10; reg7=reg17*reg35;
   reg28=reg37+reg28; reg32=reg17*reg32; reg35=reg29*reg35; reg5=reg26+reg5; reg32=reg2+reg32;
   reg6=reg6*reg31; reg12=reg12*reg33; reg30=reg40/reg30; reg31=reg16*reg31; reg33=reg22*reg33;
   reg38=reg39+reg38; reg2=reg25*reg34; reg1=reg8+reg1; reg10=reg9+reg10; reg34=reg29*reg34;
   reg24=reg27+reg24; reg35=reg20+reg35; reg7=reg28+reg7; reg0=reg25*reg0; reg12=reg6+reg12;
   reg23=reg23*reg30; reg38=reg32*reg38; reg1=reg4*reg1; reg7=reg3*reg7; reg0=reg10+reg0;
   reg35=reg32*reg35; reg34=reg24+reg34; reg33=reg31+reg33; reg2=reg5+reg2; reg30=reg25*reg30;
   reg23=reg12+reg23; reg30=reg33+reg30; reg1=reg38-reg1; reg2=reg13*reg2; reg7=reg35-reg7;
   reg34=reg0*reg34; reg3=0.011164549684630114537*reg1; reg4=0.1555021169820365473*reg1; reg5=0.1555021169820365473*reg7; reg6=0.04166666666666666908*reg7;
   reg2=reg34-reg2; reg30=reg14*reg30; reg8=0.04166666666666666908*reg1; reg23=reg0*reg23; reg0=0.011164549684630114537*reg7;
   reg30=reg23-reg30; reg4=reg6+reg4; reg9=0.04166666666666666908*reg2; reg3=reg6+reg3; reg6=0.1555021169820365473*reg2;
   reg10=0.011164549684630114537*reg2; reg0=reg0+reg8; reg5=reg8+reg5; reg6=reg3+reg6; reg3=0.011164549684630114537*reg30;
   reg5=reg9+reg5; reg10=reg4+reg10; reg4=0.04166666666666666908*reg30; reg8=0.1555021169820365473*reg30; reg9=reg0+reg9;
   Ne(0,9)+=reg4+reg10; Ne(1,10)+=reg4+reg10; Ne(2,11)+=reg4+reg10; Ne(0,6)+=reg5+reg3; Ne(1,7)+=reg5+reg3;
   Ne(2,8)+=reg5+reg3; Ne(0,3)+=reg6+reg4; Ne(1,4)+=reg6+reg4; Ne(2,5)+=reg6+reg4; Ne(0,0)+=reg9+reg8;
   Ne(1,1)+=reg9+reg8; Ne(2,2)+=reg9+reg8;

}
 /// pour les Bar_3
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Bar_3,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=2.1547005383792514621*e.pos(1)[1]; double reg1=2.1547005383792514621*e.pos(0)[1]; double reg2=0.15470053837925146212*e.pos(1)[1]; double reg3=0.15470053837925146212*e.pos(0)[0]; double reg4=0.15470053837925146212*e.pos(0)[1];
   double reg5=2.1547005383792514621*e.pos(0)[0]; double reg6=0.15470053837925146212*e.pos(1)[0]; double reg7=2.1547005383792514621*e.pos(1)[0]; reg0=reg0+reg4; reg7=reg7+reg3;
   reg5=reg6+reg5; double reg8=2.3094010767585029242*e.pos(2)[1]; reg1=reg2+reg1; double reg9=2.3094010767585029242*e.pos(2)[0]; reg7=reg7-reg9;
   reg0=reg0-reg8; double reg10=reg8-reg1; double reg11=reg9-reg5; double reg12=pow(reg11,2); double reg13=pow(reg7,2);
   double reg14=pow(reg10,2); double reg15=pow(reg0,2); reg15=reg13+reg15; reg14=reg12+reg14; reg15=pow(reg15,0.5);
   reg14=pow(reg14,0.5); reg12=reg7/reg15; reg13=reg10/reg14; reg15=reg0/reg15; reg14=reg11/reg14;
   reg12=reg7*reg12; reg15=reg0*reg15; reg13=reg10*reg13; reg14=reg11*reg14; reg13=reg14+reg13;
   reg15=reg12+reg15; reg0=0.22767090063073975644*reg13; reg7=0.061004233964073109089*reg15; reg10=0.061004233964073109089*reg13; reg11=0.22767090063073975644*reg15;
   reg12=0.33333333333333335264*reg13; reg14=0.33333333333333335264*reg15; Ne(0,0)+=reg0-reg7; Ne(1,1)+=reg0-reg7; Ne(0,2)+=reg11-reg10;
   Ne(1,3)+=reg11-reg10; Ne(0,4)+=reg12+reg14; Ne(1,5)+=reg12+reg14;

}
 /// pour les Triangle_6
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Triangle_6,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.56758792732771999991*e.pos(1)[0]; double reg3=0.78379396366385999995*e.pos(0)[0]; double reg4=0.78379396366385999995*e.pos(1)[1];
   double reg5=0.56758792732771999991*e.pos(0)[1]; double reg6=0.78379396366385999995*e.pos(0)[1]; double reg7=0.56758792732771999991*e.pos(1)[1]; double reg8=1.3513818909915799999*e.pos(3)[1]; double reg9=reg4+reg5;
   double reg10=reg0+reg1; double reg11=reg2+reg3; double reg12=0.78379396366385999995*e.pos(1)[2]; double reg13=1.3513818909915799999*e.pos(3)[0]; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg7+reg6; double reg16=0.56758792732771999991*e.pos(1)[2]; double reg17=0.78379396366385999995*e.pos(0)[2]; double reg18=1.3513818909915799999*e.pos(3)[2]; double reg19=1.78379396366386*e.pos(4)[1];
   double reg20=reg13-reg11; double reg21=2.2673902919218320001*e.pos(0)[0]; double reg22=reg12+reg14; reg9=reg9-reg8; double reg23=2.2673902919218320001*e.pos(0)[1];
   double reg24=0.63369514596091600003*e.pos(1)[1]; double reg25=reg8-reg15; double reg26=0.63369514596091600003*e.pos(1)[0]; double reg27=1.78379396366386*e.pos(4)[0]; double reg28=reg16+reg17;
   reg10=reg10-reg13; double reg29=2.9010854378827480001*e.pos(3)[1]; reg0=reg0-reg3; double reg30=0.63369514596091600003*e.pos(0)[0]; double reg31=0.63369514596091600003*e.pos(0)[1];
   reg20=reg27+reg20; double reg32=reg24+reg23; double reg33=2.2673902919218320001*e.pos(1)[1]; double reg34=2.2673902919218320001*e.pos(0)[2]; double reg35=0.63369514596091600003*e.pos(1)[2];
   reg4=reg4-reg6; double reg36=0.43241207267228000009*e.pos(4)[1]; double reg37=2.2673902919218320001*e.pos(1)[0]; double reg38=reg26+reg21; double reg39=reg18-reg28;
   double reg40=0.43241207267228000009*e.pos(4)[0]; reg10=reg27+reg10; double reg41=1.78379396366386*e.pos(5)[0]; double reg42=2.9010854378827480001*e.pos(3)[0]; reg9=reg9+reg19;
   double reg43=1.78379396366386*e.pos(5)[1]; reg25=reg19+reg25; reg22=reg22-reg18; double reg44=1.78379396366386*e.pos(4)[2]; double reg45=reg42-reg38;
   double reg46=reg29-reg32; double reg47=0.366304854039084*e.pos(4)[1]; reg12=reg12-reg17; double reg48=0.43241207267228000009*e.pos(5)[1]; reg4=reg36+reg4;
   double reg49=0.366304854039084*e.pos(4)[0]; double reg50=reg35+reg34; reg39=reg44+reg39; reg10=reg10-reg41; reg25=reg25-reg43;
   reg43=reg9-reg43; reg22=reg22+reg44; reg9=1.78379396366386*e.pos(5)[2]; reg41=reg20-reg41; reg20=0.63369514596091600003*e.pos(0)[2];
   double reg51=2.2673902919218320001*e.pos(1)[2]; reg33=reg33+reg31; reg37=reg37+reg30; reg0=reg40+reg0; double reg52=0.43241207267228000009*e.pos(4)[2];
   double reg53=0.43241207267228000009*e.pos(5)[0]; double reg54=2.9010854378827480001*e.pos(3)[2]; reg12=reg52+reg12; reg37=reg37-reg42; double reg55=0.43241207267228000009*e.pos(5)[2];
   reg48=reg4-reg48; reg53=reg0-reg53; reg0=pow(reg41,2); reg4=pow(reg25,2); reg39=reg39-reg9;
   double reg56=0.366304854039084*e.pos(4)[2]; double reg57=reg30-reg26; reg9=reg22-reg9; reg51=reg51+reg20; reg22=pow(reg43,2);
   reg33=reg33-reg29; double reg58=0.78379396366385999995*e.pos(2)[0]; double reg59=pow(reg10,2); double reg60=2.710505431213761085e-20*e.pos(3)[0]; double reg61=reg31-reg24;
   double reg62=0.78379396366385999995*e.pos(2)[1]; double reg63=reg54-reg50; double reg64=0.366304854039084*e.pos(5)[1]; double reg65=2.710505431213761085e-20*e.pos(3)[1]; reg46=reg46+reg47;
   double reg66=0.366304854039084*e.pos(5)[0]; reg45=reg49+reg45; double reg67=2.710505431213761085e-20*e.pos(3)[2]; reg57=reg57-reg60; reg61=reg61-reg65;
   double reg68=reg20-reg35; reg55=reg12-reg55; reg12=3.2673902919218320001*e.pos(4)[1]; double reg69=3.2673902919218320001*e.pos(4)[0]; reg51=reg51-reg54;
   reg33=reg47+reg33; reg37=reg49+reg37; reg46=reg46-reg64; double reg70=0.366304854039084*e.pos(5)[2]; reg63=reg63+reg56;
   double reg71=pow(reg53,2); reg45=reg45-reg66; double reg72=pow(reg48,2); double reg73=0.56758792732771999991*e.pos(2)[1]; double reg74=0.56758792732771999991*e.pos(2)[0];
   double reg75=0.78379396366385999995*e.pos(2)[2]; double reg76=reg1+reg58; double reg77=1.78379396366386*e.pos(3)[0]; double reg78=pow(reg39,2); double reg79=reg5+reg62;
   reg4=reg0+reg4; reg22=reg59+reg22; reg0=1.78379396366386*e.pos(3)[1]; reg59=pow(reg9,2); reg57=reg69+reg57;
   double reg80=1.78379396366386*e.pos(3)[2]; double reg81=3.2673902919218320001*e.pos(5)[1]; double reg82=pow(reg45,2); double reg83=pow(reg55,2); double reg84=3.2673902919218320001*e.pos(4)[2];
   reg51=reg56+reg51; double reg85=reg14+reg75; reg61=reg12+reg61; reg59=reg22+reg59; reg76=reg76-reg77;
   reg22=reg3+reg74; double reg86=0.56758792732771999991*e.pos(2)[2]; double reg87=reg6+reg73; double reg88=3.2673902919218320001*e.pos(5)[0]; reg63=reg63-reg70;
   reg78=reg4+reg78; reg4=0.63369514596091600003*e.pos(2)[1]; reg64=reg33-reg64; reg33=0.63369514596091600003*e.pos(2)[0]; reg68=reg68-reg67;
   reg3=reg58-reg3; reg58=0.43241207267228000009*e.pos(3)[0]; reg79=reg79-reg0; reg66=reg37-reg66; reg6=reg62-reg6;
   reg37=0.43241207267228000009*e.pos(3)[1]; reg72=reg71+reg72; reg62=pow(reg46,2); reg70=reg51-reg70; reg83=reg72+reg83;
   reg68=reg84+reg68; reg88=reg57-reg88; reg51=3.2673902919218320001*e.pos(5)[2]; reg77=reg77+reg22; reg81=reg61-reg81;
   reg0=reg0+reg87; reg57=reg17+reg86; reg61=0.63369514596091600003*e.pos(2)[2]; reg76=reg27+reg76; reg23=reg23+reg4;
   reg59=pow(reg59,0.5); reg85=reg85-reg80; reg62=reg82+reg62; reg71=1.3513818909915799999*e.pos(5)[1]; reg79=reg19+reg79;
   reg72=pow(reg63,2); reg82=0.366304854039084*e.pos(3)[1]; reg21=reg21+reg33; double reg89=0.366304854039084*e.pos(3)[0]; reg78=pow(reg78,0.5);
   double reg90=0.43241207267228000009*e.pos(3)[2]; double reg91=pow(reg66,2); double reg92=pow(reg64,2); reg17=reg75-reg17; reg75=1.3513818909915799999*e.pos(5)[0];
   reg58=reg3-reg58; reg37=reg6-reg37; reg3=pow(reg70,2); reg90=reg17-reg90; reg6=reg23+reg82;
   reg17=reg30-reg33; double reg93=3.2673902919218320001*e.pos(3)[0]; reg80=reg80+reg57; double reg94=reg31-reg4; double reg95=reg10/reg59;
   double reg96=3.2673902919218320001*e.pos(3)[1]; reg51=reg68-reg51; reg68=1.3513818909915799999*e.pos(5)[2]; reg85=reg44+reg85; reg36=reg37+reg36;
   reg79=reg79-reg71; reg58=reg40+reg58; reg37=reg43/reg59; reg76=reg76-reg75; reg92=reg91+reg92;
   reg72=reg62+reg72; reg40=reg25/reg78; reg62=reg41/reg78; reg91=reg21+reg89; double reg97=0.366304854039084*e.pos(3)[2];
   double reg98=2.2673902919218320001*e.pos(2)[1]; reg83=pow(reg83,0.5); reg19=reg19-reg0; double reg99=pow(reg81,2); reg27=reg27-reg77;
   double reg100=pow(reg88,2); double reg101=2.2673902919218320001*e.pos(2)[0]; reg34=reg34+reg61; double reg102=reg48/reg83; double reg103=reg37*reg79;
   reg96=reg94-reg96; reg85=reg85-reg68; reg94=2.9010854378827480001*e.pos(5)[1]; double reg104=reg20-reg61; double reg105=3.2673902919218320001*e.pos(3)[2];
   reg52=reg90+reg52; reg90=reg40*reg36; double reg106=reg47-reg6; double reg107=reg53/reg83; double reg108=reg62*reg58;
   reg101=reg30+reg101; double reg109=reg95*reg76; double reg110=reg34+reg97; reg78=reg39/reg78; reg99=reg100+reg99;
   reg72=pow(reg72,0.5); reg59=reg9/reg59; reg100=2.9010854378827480001*e.pos(5)[0]; double reg111=reg49-reg91; double reg112=pow(reg51,2);
   reg19=reg71+reg19; double reg113=2.2673902919218320001*e.pos(2)[2]; reg93=reg17-reg93; reg27=reg75+reg27; reg3=reg92+reg3;
   reg98=reg31+reg98; reg44=reg44-reg80; reg17=2.9010854378827480001*e.pos(5)[2]; reg92=reg45/reg72; double reg114=reg59*reg85;
   reg83=reg55/reg83; double reg115=reg107*reg27; reg90=reg108+reg90; reg44=reg68+reg44; reg103=reg109+reg103;
   reg113=reg20+reg113; reg82=reg98-reg82; reg93=reg69+reg93; reg89=reg101-reg89; reg69=reg46/reg72;
   reg98=reg56-reg110; reg111=reg100+reg111; reg101=reg102*reg19; reg112=reg99+reg112; reg3=pow(reg3,0.5);
   reg105=reg104-reg105; reg99=reg78*reg52; reg106=reg106+reg94; reg12=reg96+reg12; reg96=5.42101086242752217e-20*e.pos(5)[0];
   reg104=5.42101086242752217e-20*e.pos(5)[1]; reg101=reg115+reg101; reg112=pow(reg112,0.5); reg72=reg63/reg72; reg108=reg83*reg44;
   reg82=reg47+reg82; reg12=reg12-reg104; reg93=reg93-reg96; reg47=reg92*reg111; reg109=reg66/reg3;
   reg115=reg64/reg3; double reg116=reg69*reg106; reg89=reg49+reg89; reg97=reg113-reg97; reg114=reg103+reg114;
   reg49=5.42101086242752217e-20*e.pos(5)[2]; reg84=reg105+reg84; reg99=reg90+reg99; reg98=reg98+reg17; reg116=reg47+reg116;
   reg84=reg84-reg49; reg3=reg70/reg3; reg47=reg88/reg112; reg90=reg81/reg112; reg103=reg62*reg99;
   reg89=reg89-reg100; reg105=reg37*reg114; reg113=reg115*reg12; double reg117=reg72*reg98; reg82=reg82-reg94;
   double reg118=reg40*reg99; double reg119=reg95*reg114; reg97=reg56+reg97; reg108=reg101+reg108; reg56=reg109*reg93;
   reg117=reg116+reg117; reg113=reg56+reg113; reg119=reg76-reg119; reg105=reg79-reg105; reg56=reg59*reg114;
   reg101=reg3*reg84; reg116=reg102*reg108; double reg120=reg107*reg108; reg118=reg36-reg118; double reg121=reg78*reg99;
   reg112=reg51/reg112; double reg122=reg47*reg89; double reg123=reg90*reg82; reg103=reg58-reg103; reg97=reg97-reg17;
   double reg124=reg69*reg117; double reg125=reg83*reg108; double reg126=reg92*reg117; double reg127=pow(reg119,2); double reg128=pow(reg105,2);
   reg56=reg85-reg56; reg101=reg113+reg101; reg113=pow(reg118,2); reg116=reg19-reg116; reg120=reg27-reg120;
   reg121=reg52-reg121; reg123=reg122+reg123; reg122=pow(reg103,2); double reg129=reg112*reg97; double reg130=reg109*reg101;
   reg129=reg123+reg129; reg125=reg44-reg125; reg123=pow(reg116,2); double reg131=pow(reg120,2); double reg132=reg115*reg101;
   reg122=reg113+reg122; reg113=pow(reg121,2); double reg133=pow(reg56,2); reg128=reg127+reg128; reg126=reg111-reg126;
   reg124=reg106-reg124; reg127=reg72*reg117; reg127=reg98-reg127; reg113=reg122+reg113; reg122=reg47*reg129;
   double reg134=reg90*reg129; reg123=reg131+reg123; reg131=pow(reg125,2); double reg135=reg3*reg101; reg132=reg12-reg132;
   reg130=reg93-reg130; reg133=reg128+reg133; reg128=pow(reg126,2); double reg136=pow(reg124,2); reg135=reg84-reg135;
   double reg137=pow(reg132,2); double reg138=pow(reg130,2); double reg139=reg112*reg129; reg122=reg89-reg122; reg134=reg82-reg134;
   reg131=reg123+reg131; reg133=pow(reg133,0.5); reg136=reg128+reg136; reg123=pow(reg127,2); reg113=pow(reg113,0.5);
   reg128=pow(reg135,2); reg131=pow(reg131,0.5); reg118=reg118/reg113; reg119=reg119/reg133; reg105=reg105/reg133;
   reg103=reg103/reg113; reg123=reg136+reg123; reg137=reg138+reg137; reg136=pow(reg134,2); reg139=reg97-reg139;
   reg138=pow(reg122,2); double reg140=pow(reg139,2); double reg141=reg41*reg103; reg37=reg43*reg37; reg76=reg76*reg119;
   reg36=reg36*reg118; reg103=reg58*reg103; reg95=reg10*reg95; reg136=reg138+reg136; reg113=reg121/reg113;
   reg120=reg120/reg131; reg62=reg41*reg62; reg40=reg25*reg40; reg123=pow(reg123,0.5); reg43=reg43*reg105;
   reg116=reg116/reg131; reg119=reg10*reg119; reg118=reg25*reg118; reg128=reg137+reg128; reg133=reg56/reg133;
   reg105=reg79*reg105; reg102=reg48*reg102; reg27=reg27*reg120; reg128=pow(reg128,0.5); reg19=reg19*reg116;
   reg37=reg95+reg37; reg131=reg125/reg131; reg120=reg53*reg120; reg52=reg52*reg113; reg113=reg39*reg113;
   reg118=reg141+reg118; reg85=reg85*reg133; reg59=reg9*reg59; reg140=reg136+reg140; reg133=reg9*reg133;
   reg43=reg119+reg43; reg40=reg62+reg40; reg78=reg39*reg78; reg116=reg48*reg116; reg105=reg76+reg105;
   reg124=reg124/reg123; reg107=reg53*reg107; reg126=reg126/reg123; reg36=reg103+reg36; reg140=pow(reg140,0.5);
   reg85=reg105+reg85; reg113=reg118+reg113; reg19=reg27+reg19; reg102=reg107+reg102; reg133=reg43+reg133;
   reg83=reg55*reg83; reg44=reg44*reg131; reg106=reg106*reg124; reg111=reg111*reg126; reg123=reg127/reg123;
   reg69=reg46*reg69; reg92=reg45*reg92; reg78=reg40+reg78; reg131=reg55*reg131; reg116=reg120+reg116;
   reg126=reg45*reg126; reg130=reg130/reg128; reg124=reg46*reg124; reg132=reg132/reg128; reg52=reg36+reg52;
   reg59=reg37+reg59; reg115=reg64*reg115; reg69=reg92+reg69; reg109=reg66*reg109; reg133=reg114*reg133;
   reg113=reg99*reg113; reg85=reg59*reg85; reg124=reg126+reg124; reg98=reg98*reg123; reg123=reg63*reg123;
   reg122=reg122/reg140; reg66=reg66*reg130; reg64=reg64*reg132; reg128=reg135/reg128; reg132=reg12*reg132;
   reg134=reg134/reg140; reg130=reg93*reg130; reg106=reg111+reg106; reg44=reg19+reg44; reg131=reg116+reg131;
   reg83=reg102+reg83; reg52=reg78*reg52; reg72=reg63*reg72; reg64=reg66+reg64; reg113=reg52-reg113;
   reg132=reg130+reg132; reg3=reg70*reg3; reg9=reg88*reg122; reg115=reg109+reg115; reg123=reg124+reg123;
   reg98=reg106+reg98; reg72=reg69+reg72; reg131=reg108*reg131; reg44=reg83*reg44; reg84=reg84*reg128;
   reg140=reg139/reg140; reg128=reg70*reg128; reg10=reg81*reg134; reg134=reg82*reg134; reg122=reg89*reg122;
   reg133=reg85-reg133; reg47=reg88*reg47; reg90=reg81*reg90; reg10=reg9+reg10; reg123=reg117*reg123;
   reg9=reg51*reg140; reg12=0.005384432036113586778*reg133; reg19=0.088847818743090689935*reg113; reg25=0.009463616120767210603*reg113; reg27=0.021537728144454347112*reg133;
   reg36=0.021537728144454347112*reg113; reg140=reg97*reg140; reg134=reg122+reg134; reg37=0.088847818743090689935*reg133; reg112=reg51*reg112;
   reg90=reg47+reg90; reg39=0.009463616120767210603*reg133; reg131=reg44-reg131; reg128=reg64+reg128; reg40=0.005384432036113586778*reg113;
   reg84=reg132+reg84; reg3=reg115+reg3; reg98=reg72*reg98; reg9=reg10+reg9; reg19=reg27+reg19;
   reg123=reg98-reg123; reg25=reg12+reg25; reg37=reg36+reg37; reg112=reg90+reg112; reg10=0.021537728144454347112*reg131;
   reg12=reg40+reg12; reg41=0.009463616120767210603*reg131; reg43=0.005384432036113586778*reg131; reg36=reg27+reg36; reg40=reg39+reg40;
   reg140=reg134+reg140; reg84=reg3*reg84; reg128=reg101*reg128; reg3=0.088847818743090689935*reg131; reg40=reg40+reg43;
   reg41=reg12+reg41; reg3=reg36+reg3; reg12=0.016449618187943419918*reg123; reg37=reg37+reg10; reg128=reg84-reg128;
   reg140=reg112*reg140; reg27=0.028457289286966203713*reg123; reg9=reg129*reg9; reg19=reg10+reg19; reg25=reg43+reg25;
   reg10=0.0041124045469858549794*reg123; reg36=0.0018441552587796664112*reg123; reg19=reg12+reg19; reg39=0.016449618187943419918*reg128; reg36=reg37+reg36;
   reg37=0.0018441552587796664109*reg128; reg43=0.0041124045469858549794*reg128; reg44=0.016449618187943419916*reg128; reg45=0.004112404546985854979*reg128; reg40=reg27-reg40;
   reg9=reg140-reg9; reg27=0.028457289286966203713*reg128; reg25=reg25+reg10; reg41=reg10+reg41; reg12=reg3+reg12;
   reg43=reg41+reg43; reg3=0.028457289286966203713*reg9; reg25=reg27-reg25; reg37=reg19+reg37; reg45=reg40-reg45;
   reg39=reg36+reg39; reg10=0.016449618187943419918*reg9; reg44=reg12+reg44; reg12=0.0041124045469858549794*reg9; reg19=0.0018441552587796664111*reg9;
   Ne(0,15)+=reg10+reg37; Ne(1,16)+=reg10+reg37; Ne(2,17)+=reg10+reg37; Ne(0,0)+=reg45-reg12; Ne(1,1)+=reg45-reg12;
   Ne(2,2)+=reg45-reg12; Ne(0,12)+=reg39+reg10; Ne(1,13)+=reg39+reg10; Ne(2,14)+=reg39+reg10; Ne(0,9)+=reg44+reg19;
   Ne(1,10)+=reg44+reg19; Ne(2,11)+=reg44+reg19; Ne(0,6)+=reg3-reg43; Ne(1,7)+=reg3-reg43; Ne(2,8)+=reg3-reg43;
   Ne(0,3)+=reg25-reg12; Ne(1,4)+=reg25-reg12; Ne(2,5)+=reg25-reg12;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(1)[0]; double reg1=0.12200846792814621817*e.pos(1)[0]; double reg2=0.36602540378443865451*e.pos(0)[0]; double reg3=0.12200846792814621817*e.pos(0)[1]; double reg4=0.36602540378443865451*e.pos(1)[1];
   double reg5=0.12200846792814621817*e.pos(0)[0]; double reg6=0.12200846792814621817*e.pos(1)[1]; double reg7=0.36602540378443865451*e.pos(0)[1]; double reg8=0.12200846792814621817*e.pos(1)[2]; double reg9=0.36602540378443865451*e.pos(0)[2];
   double reg10=reg1+reg2; double reg11=1.3660254037844385386*e.pos(2)[0]; double reg12=0.45534180126147951289*e.pos(2)[1]; double reg13=reg4+reg3; double reg14=1.3660254037844385386*e.pos(2)[1];
   double reg15=reg6+reg7; double reg16=0.45534180126147951289*e.pos(2)[0]; double reg17=reg0+reg5; double reg18=0.36602540378443865451*e.pos(1)[2]; double reg19=0.12200846792814621817*e.pos(0)[2];
   double reg20=1.3660254037844385386*e.pos(1)[0]; double reg21=0.45534180126147951289*e.pos(0)[0]; double reg22=0.45534180126147951289*e.pos(0)[1]; double reg23=1.3660254037844385386*e.pos(1)[1]; double reg24=0.45534180126147951289*e.pos(1)[0];
   double reg25=1.3660254037844385386*e.pos(0)[0]; double reg26=1.3660254037844385386*e.pos(2)[2]; double reg27=reg21+reg20; double reg28=reg22+reg23; double reg29=1.3660254037844385386*e.pos(1)[2];
   double reg30=reg8+reg9; double reg31=1.3660254037844385386*e.pos(0)[1]; double reg32=1.3660254037844385386*e.pos(3)[1]; double reg33=0.45534180126147951289*e.pos(1)[1]; double reg34=0.45534180126147951289*e.pos(0)[2];
   double reg35=reg17+reg16; double reg36=1.3660254037844385386*e.pos(3)[0]; double reg37=0.12200846792814621817*e.pos(2)[1]; double reg38=reg13+reg12; double reg39=0.12200846792814621817*e.pos(2)[0];
   double reg40=0.45534180126147951289*e.pos(2)[2]; double reg41=reg18+reg19; double reg42=0.45534180126147951289*e.pos(3)[0]; reg10=reg10+reg11; double reg43=0.45534180126147951289*e.pos(3)[1];
   reg15=reg15+reg14; double reg44=reg35+reg36; double reg45=reg34+reg29; double reg46=1.3660254037844385386*e.pos(3)[2]; double reg47=reg41+reg40;
   double reg48=0.48803387171258487271*e.pos(4)[1]; double reg49=reg24+reg25; double reg50=0.36602540378443865451*e.pos(2)[0]; double reg51=0.36602540378443865451*e.pos(2)[1]; reg28=reg37+reg28;
   double reg52=0.48803387171258487271*e.pos(4)[0]; double reg53=0.12200846792814621817*e.pos(2)[2]; reg10=reg10+reg42; double reg54=0.45534180126147951289*e.pos(3)[2]; reg30=reg30+reg26;
   double reg55=reg38+reg32; reg15=reg15+reg43; reg27=reg39+reg27; double reg56=0.36602540378443865451*e.pos(3)[0]; double reg57=0.36602540378443865451*e.pos(3)[1];
   double reg58=1.3660254037844385386*e.pos(0)[2]; double reg59=0.45534180126147951289*e.pos(1)[2]; double reg60=reg33+reg31; reg10=reg10-reg52; double reg61=reg52-reg44;
   double reg62=reg60+reg51; double reg63=0.12200846792814621817*e.pos(3)[1]; double reg64=0.66666666666666670528*e.pos(5)[0]; double reg65=1.8213672050459180516*e.pos(4)[0]; reg27=reg27+reg56;
   double reg66=0.36602540378443865451*e.pos(3)[2]; reg15=reg15-reg48; double reg67=0.66666666666666670528*e.pos(5)[1]; double reg68=1.8213672050459180516*e.pos(4)[1]; reg28=reg28+reg57;
   double reg69=reg59+reg58; reg30=reg30+reg54; double reg70=0.48803387171258487271*e.pos(4)[2]; reg45=reg53+reg45; double reg71=0.12200846792814621817*e.pos(3)[0];
   double reg72=reg47+reg46; double reg73=reg49+reg50; double reg74=0.36602540378443865451*e.pos(2)[2]; double reg75=reg48-reg55; reg45=reg45+reg66;
   double reg76=0.12200846792814621817*e.pos(3)[2]; double reg77=1.8213672050459180516*e.pos(4)[2]; reg27=reg27-reg65; reg28=reg28-reg68; double reg78=reg69+reg74;
   reg15=reg15+reg67; double reg79=1.8213672050459180516*e.pos(6)[1]; double reg80=reg63+reg62; reg30=reg30-reg70; double reg81=0.66666666666666670528*e.pos(5)[2];
   double reg82=reg70-reg72; double reg83=reg71+reg73; reg75=reg67+reg75; reg61=reg64+reg61; reg10=reg10+reg64;
   double reg84=1.8213672050459180516*e.pos(6)[0]; double reg85=reg0+reg21; double reg86=0.48803387171258487271*e.pos(6)[0]; reg75=reg79+reg75; reg82=reg81+reg82;
   reg61=reg84+reg61; reg27=reg64+reg27; double reg87=reg68-reg80; reg10=reg10-reg84; reg45=reg45-reg77;
   reg28=reg67+reg28; reg15=reg15-reg79; double reg88=0.66666666666666670528*e.pos(7)[1]; reg30=reg30+reg81; double reg89=1.8213672050459180516*e.pos(6)[2];
   double reg90=reg65-reg83; double reg91=reg2+reg24; double reg92=0.66666666666666670528*e.pos(7)[0]; double reg93=reg7+reg33; double reg94=reg4+reg22;
   double reg95=reg76+reg78; double reg96=0.48803387171258487271*e.pos(6)[1]; double reg97=reg18+reg34; reg45=reg81+reg45; reg23=reg3+reg23;
   reg61=reg61-reg92; reg11=reg11+reg91; reg10=reg10-reg92; reg75=reg75-reg88; reg90=reg64+reg90;
   reg27=reg27-reg86; reg3=0.66666666666666670528*e.pos(7)[2]; reg30=reg30-reg89; reg82=reg89+reg82; reg15=reg15-reg88;
   reg87=reg67+reg87; reg28=reg28-reg96; reg64=0.48803387171258487271*e.pos(6)[2]; reg14=reg14+reg93; reg20=reg5+reg20;
   reg5=reg9+reg59; reg39=reg39+reg85; reg37=reg37+reg94; reg67=reg77-reg95; reg53=reg53+reg97;
   reg30=reg30-reg3; reg26=reg26+reg5; reg28=reg28-reg88; reg82=reg82-reg3; double reg98=pow(reg10,2);
   reg45=reg45-reg64; double reg99=pow(reg15,2); reg31=reg6+reg31; reg87=reg96+reg87; reg27=reg27-reg92;
   reg25=reg1+reg25; reg20=reg16+reg20; reg63=reg14+reg63; reg1=pow(reg61,2); reg23=reg12+reg23;
   reg67=reg81+reg67; reg6=0.66666666666666670528*e.pos(4)[1]; reg14=0.66666666666666670528*e.pos(4)[0]; reg71=reg11+reg71; reg37=reg32+reg37;
   reg39=reg36+reg39; reg11=pow(reg75,2); reg29=reg19+reg29; reg90=reg86+reg90; reg39=reg39-reg14;
   reg19=0.48803387171258487271*e.pos(5)[0]; reg67=reg64+reg67; reg32=0.48803387171258487271*e.pos(5)[1]; reg53=reg46+reg53; reg37=reg37-reg6;
   reg36=0.66666666666666670528*e.pos(4)[2]; reg71=reg71-reg14; reg46=1.8213672050459180516*e.pos(5)[0]; reg29=reg40+reg29; reg81=reg57+reg23;
   reg63=reg63-reg6; double reg100=1.8213672050459180516*e.pos(5)[1]; double reg101=reg56+reg20; reg92=reg90-reg92; reg90=pow(reg30,2);
   reg76=reg26+reg76; reg99=reg98+reg99; reg45=reg45-reg3; reg26=pow(reg28,2); reg58=reg8+reg58;
   reg8=pow(reg27,2); reg31=reg51+reg31; reg11=reg1+reg11; reg1=pow(reg82,2); reg25=reg50+reg25;
   reg88=reg87-reg88; reg87=pow(reg88,2); reg98=reg66+reg29; reg53=reg53-reg36; double reg102=0.48803387171258487271*e.pos(5)[2];
   reg3=reg67-reg3; reg67=reg6+reg81; double reg103=reg14+reg101; double reg104=pow(reg92,2); double reg105=pow(reg45,2);
   reg26=reg8+reg26; reg76=reg76-reg36; reg63=reg63-reg100; reg8=0.66666666666666670528*e.pos(6)[1]; reg58=reg74+reg58;
   double reg106=reg43+reg31; reg71=reg71-reg46; reg1=reg11+reg1; reg11=reg42+reg25; reg39=reg39-reg19;
   reg37=reg37-reg32; double reg107=0.66666666666666670528*e.pos(6)[0]; double reg108=1.8213672050459180516*e.pos(5)[2]; reg90=reg99+reg90; reg39=reg107+reg39;
   reg71=reg107+reg71; reg99=0.48803387171258487271*e.pos(7)[0]; reg14=reg14+reg11; reg105=reg26+reg105; reg26=reg100-reg67;
   reg1=pow(reg1,0.5); reg90=pow(reg90,0.5); reg6=reg6+reg106; reg87=reg104+reg87; reg104=0.48803387171258487271*e.pos(7)[1];
   double reg109=1.8213672050459180516*e.pos(7)[0]; double reg110=pow(reg3,2); double reg111=0.66666666666666670528*e.pos(6)[2]; reg76=reg76-reg108; double reg112=reg36+reg98;
   double reg113=reg46-reg103; reg53=reg53-reg102; reg63=reg63+reg8; double reg114=1.8213672050459180516*e.pos(7)[1]; reg37=reg8+reg37;
   double reg115=reg54+reg58; double reg116=reg32-reg6; reg63=reg63-reg104; double reg117=1.8213672050459180516*e.pos(7)[2]; reg36=reg36+reg115;
   reg113=reg107+reg113; reg71=reg71-reg99; reg105=pow(reg105,0.5); reg53=reg111+reg53; double reg118=reg19-reg14;
   double reg119=reg61/reg1; double reg120=reg75/reg1; reg26=reg8+reg26; double reg121=reg10/reg90; reg39=reg39-reg109;
   double reg122=reg15/reg90; reg110=reg87+reg110; reg87=0.48803387171258487271*e.pos(7)[2]; reg37=reg37-reg114; reg76=reg76+reg111;
   double reg123=reg108-reg112; double reg124=reg122*reg63; double reg125=reg28/reg105; double reg126=reg27/reg105; reg76=reg76-reg87;
   double reg127=reg121*reg71; double reg128=reg120*reg37; reg53=reg53-reg117; reg90=reg30/reg90; reg123=reg111+reg123;
   reg110=pow(reg110,0.5); double reg129=reg119*reg39; reg26=reg104+reg26; reg1=reg82/reg1; reg118=reg107+reg118;
   reg116=reg8+reg116; reg113=reg99+reg113; reg8=reg102-reg36; reg124=reg127+reg124; reg105=reg45/reg105;
   reg107=reg125*reg26; reg123=reg87+reg123; reg127=reg126*reg113; reg128=reg129+reg128; reg129=reg90*reg76;
   double reg130=reg1*reg53; double reg131=reg92/reg110; double reg132=reg88/reg110; reg118=reg109+reg118; reg116=reg114+reg116;
   reg8=reg111+reg8; reg107=reg127+reg107; reg129=reg124+reg129; reg8=reg117+reg8; reg111=reg105*reg123;
   reg110=reg3/reg110; reg124=reg132*reg116; reg127=reg131*reg118; reg130=reg128+reg130; reg128=reg121*reg129;
   double reg133=reg120*reg130; reg111=reg107+reg111; reg107=reg119*reg130; double reg134=reg122*reg129; reg124=reg127+reg124;
   reg127=reg110*reg8; double reg135=reg126*reg111; double reg136=reg125*reg111; reg133=reg37-reg133; reg128=reg71-reg128;
   reg127=reg124+reg127; reg124=reg1*reg130; reg107=reg39-reg107; reg134=reg63-reg134; double reg137=reg90*reg129;
   double reg138=reg132*reg127; reg136=reg26-reg136; double reg139=reg105*reg111; double reg140=reg131*reg127; double reg141=pow(reg107,2);
   reg135=reg113-reg135; double reg142=pow(reg133,2); double reg143=pow(reg134,2); reg124=reg53-reg124; double reg144=pow(reg128,2);
   reg137=reg76-reg137; reg140=reg118-reg140; double reg145=reg110*reg127; reg138=reg116-reg138; double reg146=pow(reg124,2);
   double reg147=pow(reg137,2); reg142=reg141+reg142; reg141=pow(reg135,2); double reg148=pow(reg136,2); reg143=reg144+reg143;
   reg139=reg123-reg139; reg147=reg143+reg147; reg143=pow(reg139,2); reg145=reg8-reg145; reg144=pow(reg138,2);
   reg148=reg141+reg148; reg141=pow(reg140,2); reg146=reg142+reg146; reg143=reg148+reg143; reg146=pow(reg146,0.5);
   reg144=reg141+reg144; reg141=pow(reg145,2); reg147=pow(reg147,0.5); reg128=reg128/reg147; reg143=pow(reg143,0.5);
   reg133=reg133/reg146; reg134=reg134/reg147; reg141=reg144+reg141; reg107=reg107/reg146; reg142=reg15*reg134;
   reg144=reg10*reg128; reg146=reg124/reg146; reg39=reg39*reg107; reg124=reg75*reg133; reg107=reg61*reg107;
   reg133=reg37*reg133; reg119=reg61*reg119; reg120=reg75*reg120; reg147=reg137/reg147; reg134=reg63*reg134;
   reg128=reg71*reg128; reg135=reg135/reg143; reg122=reg15*reg122; reg121=reg10*reg121; reg141=pow(reg141,0.5);
   reg136=reg136/reg143; reg113=reg113*reg135; reg1=reg82*reg1; reg120=reg119+reg120; reg10=reg30*reg147;
   reg142=reg144+reg142; reg147=reg76*reg147; reg134=reg128+reg134; reg90=reg30*reg90; reg122=reg121+reg122;
   reg140=reg140/reg141; reg138=reg138/reg141; reg26=reg26*reg136; reg126=reg27*reg126; reg125=reg28*reg125;
   reg143=reg139/reg143; reg133=reg39+reg133; reg53=reg53*reg146; reg135=reg27*reg135; reg136=reg28*reg136;
   reg146=reg82*reg146; reg124=reg107+reg124; reg136=reg135+reg136; reg125=reg126+reg125; reg105=reg45*reg105;
   reg15=reg92*reg140; reg1=reg120+reg1; reg90=reg122+reg90; reg141=reg145/reg141; reg45=reg45*reg143;
   reg27=reg88*reg138; reg138=reg116*reg138; reg143=reg123*reg143; reg26=reg113+reg26; reg140=reg118*reg140;
   reg132=reg88*reg132; reg131=reg92*reg131; reg147=reg134+reg147; reg146=reg124+reg146; reg53=reg133+reg53;
   reg10=reg142+reg10; reg10=reg129*reg10; reg132=reg131+reg132; reg110=reg3*reg110; reg143=reg26+reg143;
   reg138=reg140+reg138; reg45=reg136+reg45; reg8=reg8*reg141; reg53=reg1*reg53; reg147=reg90*reg147;
   reg105=reg125+reg105; reg146=reg130*reg146; reg141=reg3*reg141; reg27=reg15+reg27; reg8=reg138+reg8;
   reg141=reg27+reg141; reg45=reg111*reg45; reg110=reg132+reg110; reg146=reg53-reg146; reg10=reg147-reg10;
   reg143=reg105*reg143; reg1=0.13144585576580215187*reg10; reg3=0.024056261216234395431*reg146; reg15=0.024056261216234409915*reg10; reg26=0.024056261216234395431*reg10;
   reg27=0.035220810900864524453*reg10; reg28=0.035220810900864524453*reg146; reg30=0.024056261216234409915*reg146; reg37=0.04166666666666666908*reg10; reg141=reg127*reg141;
   reg39=0.04166666666666666908*reg146; reg53=0.13144585576580215187*reg146; reg8=reg110*reg8; reg45=reg143-reg45; reg26=reg26-reg39;
   reg61=reg28+reg1; reg63=0.13144585576580215187*reg45; reg28=reg27+reg28; reg71=0.035220810900864524453*reg45; reg1=reg1+reg53;
   reg75=0.024056261216234409915*reg45; reg3=reg3-reg37; reg76=0.024056261216234395431*reg45; reg82=0.04166666666666666908*reg45; reg39=reg15+reg39;
   reg141=reg8-reg141; reg53=reg27+reg53; reg30=reg37+reg30; reg53=reg71+reg53; reg8=0.04166666666666666908*reg141;
   reg30=reg76-reg30; reg28=reg28+reg63; reg15=0.13144585576580215187*reg141; reg71=reg1+reg71; reg1=0.024056261216234395431*reg141;
   reg75=reg3-reg75; reg61=reg63+reg61; reg3=0.035220810900864524453*reg141; reg26=reg26-reg82; reg82=reg39+reg82;
   reg27=0.024056261216234409915*reg141; Ne(0,21)+=reg15+reg53; Ne(1,22)+=reg15+reg53; Ne(2,23)+=reg15+reg53; Ne(0,18)+=reg3+reg71;
   Ne(1,19)+=reg3+reg71; Ne(2,20)+=reg3+reg71; Ne(0,15)+=reg61+reg3; Ne(1,16)+=reg61+reg3; Ne(2,17)+=reg61+reg3;
   Ne(0,12)+=reg28+reg15; Ne(1,13)+=reg28+reg15; Ne(2,14)+=reg28+reg15; Ne(0,9)+=reg75-reg8; Ne(1,10)+=reg75-reg8;
   Ne(2,11)+=reg75-reg8; Ne(0,6)+=reg26-reg27; Ne(1,7)+=reg26-reg27; Ne(2,8)+=reg26-reg27; Ne(0,3)+=reg30-reg8;
   Ne(1,4)+=reg30-reg8; Ne(2,5)+=reg30-reg8; Ne(0,0)+=reg1-reg82; Ne(1,1)+=reg1-reg82; Ne(2,2)+=reg1-reg82;

}

};
} // namespace LMT

