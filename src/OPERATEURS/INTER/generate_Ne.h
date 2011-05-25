
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
   
   double reg0=0.15470053837925146212*e.pos(1)[1]; double reg1=2.1547005383792514621*e.pos(0)[1]; double reg2=2.1547005383792514621*e.pos(1)[0]; double reg3=0.15470053837925146212*e.pos(0)[0]; double reg4=2.1547005383792514621*e.pos(0)[0];
   double reg5=0.15470053837925146212*e.pos(1)[0]; double reg6=2.1547005383792514621*e.pos(1)[1]; double reg7=0.15470053837925146212*e.pos(0)[1]; reg1=reg0+reg1; reg6=reg6+reg7;
   double reg8=2.3094010767585029242*e.pos(2)[1]; reg4=reg5+reg4; double reg9=2.3094010767585029242*e.pos(2)[0]; reg2=reg2+reg3; reg6=reg6-reg8;
   reg2=reg2-reg9; double reg10=reg8-reg1; double reg11=reg9-reg4; double reg12=pow(reg6,2); double reg13=pow(reg10,2);
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
   
   double reg0=0.56758792732771999991*e.pos(1)[1]; double reg1=0.78379396366385999995*e.pos(0)[1]; double reg2=0.56758792732771999991*e.pos(0)[1]; double reg3=0.78379396366385999995*e.pos(1)[1]; double reg4=0.78379396366385999995*e.pos(0)[0];
   double reg5=0.56758792732771999991*e.pos(0)[0]; double reg6=0.78379396366385999995*e.pos(1)[0]; double reg7=0.56758792732771999991*e.pos(1)[0]; double reg8=1.3513818909915799999*e.pos(3)[1]; double reg9=reg3+reg2;
   double reg10=0.78379396366385999995*e.pos(0)[2]; double reg11=0.56758792732771999991*e.pos(1)[2]; double reg12=reg0+reg1; double reg13=0.56758792732771999991*e.pos(0)[2]; double reg14=0.78379396366385999995*e.pos(1)[2];
   double reg15=reg7+reg4; double reg16=1.3513818909915799999*e.pos(3)[0]; double reg17=reg6+reg5; double reg18=0.63369514596091600003*e.pos(1)[0]; reg9=reg9-reg8;
   double reg19=1.78379396366386*e.pos(4)[1]; double reg20=2.2673902919218320001*e.pos(0)[0]; double reg21=reg11+reg10; double reg22=2.2673902919218320001*e.pos(0)[1]; double reg23=0.63369514596091600003*e.pos(1)[1];
   double reg24=reg14+reg13; double reg25=1.78379396366386*e.pos(4)[0]; double reg26=1.3513818909915799999*e.pos(3)[2]; double reg27=reg16-reg15; reg17=reg17-reg16;
   double reg28=reg8-reg12; double reg29=0.43241207267228000009*e.pos(4)[0]; reg28=reg19+reg28; double reg30=2.2673902919218320001*e.pos(0)[2]; double reg31=reg26-reg21;
   reg27=reg25+reg27; double reg32=reg23+reg22; double reg33=2.9010854378827480001*e.pos(3)[1]; double reg34=0.63369514596091600003*e.pos(0)[0]; double reg35=0.63369514596091600003*e.pos(0)[1];
   double reg36=2.2673902919218320001*e.pos(1)[1]; reg6=reg6-reg4; reg3=reg3-reg1; reg17=reg25+reg17; double reg37=1.78379396366386*e.pos(5)[0];
   reg24=reg24-reg26; double reg38=1.78379396366386*e.pos(4)[2]; double reg39=2.2673902919218320001*e.pos(1)[0]; reg9=reg19+reg9; double reg40=0.63369514596091600003*e.pos(1)[2];
   double reg41=0.43241207267228000009*e.pos(4)[1]; double reg42=2.9010854378827480001*e.pos(3)[0]; double reg43=1.78379396366386*e.pos(5)[1]; double reg44=reg18+reg20; double reg45=2.9010854378827480001*e.pos(3)[2];
   double reg46=0.63369514596091600003*e.pos(0)[2]; double reg47=2.2673902919218320001*e.pos(1)[2]; reg39=reg39+reg34; reg6=reg29+reg6; double reg48=reg33-reg32;
   double reg49=0.43241207267228000009*e.pos(5)[0]; reg14=reg14-reg10; reg31=reg38+reg31; reg3=reg41+reg3; double reg50=0.43241207267228000009*e.pos(5)[1];
   double reg51=reg40+reg30; double reg52=0.366304854039084*e.pos(4)[1]; double reg53=0.366304854039084*e.pos(4)[0]; reg28=reg28-reg43; reg27=reg27-reg37;
   double reg54=reg42-reg44; double reg55=0.43241207267228000009*e.pos(4)[2]; reg43=reg9-reg43; reg9=1.78379396366386*e.pos(5)[2]; reg24=reg24+reg38;
   reg37=reg17-reg37; reg36=reg36+reg35; reg17=reg34-reg18; double reg56=0.78379396366385999995*e.pos(2)[1]; reg47=reg47+reg46;
   double reg57=2.710505431213761085e-20*e.pos(3)[0]; reg31=reg31-reg9; double reg58=pow(reg27,2); double reg59=0.78379396366385999995*e.pos(2)[0]; double reg60=pow(reg28,2);
   reg9=reg24-reg9; reg24=pow(reg37,2); double reg61=0.366304854039084*e.pos(4)[2]; double reg62=reg45-reg51; double reg63=pow(reg43,2);
   double reg64=reg35-reg23; double reg65=2.710505431213761085e-20*e.pos(3)[1]; reg54=reg53+reg54; double reg66=0.366304854039084*e.pos(5)[0]; reg36=reg36-reg33;
   reg50=reg3-reg50; reg39=reg39-reg42; reg49=reg6-reg49; reg14=reg55+reg14; reg48=reg52+reg48;
   reg3=0.43241207267228000009*e.pos(5)[2]; reg6=0.366304854039084*e.pos(5)[1]; reg60=reg58+reg60; reg58=2.710505431213761085e-20*e.pos(3)[2]; double reg67=pow(reg49,2);
   double reg68=pow(reg9,2); double reg69=0.366304854039084*e.pos(5)[2]; reg62=reg62+reg61; reg3=reg14-reg3; reg14=reg46-reg40;
   reg63=reg24+reg63; reg36=reg52+reg36; reg64=reg64-reg65; reg39=reg53+reg39; reg54=reg54-reg66;
   reg24=0.78379396366385999995*e.pos(2)[2]; double reg70=3.2673902919218320001*e.pos(4)[0]; double reg71=3.2673902919218320001*e.pos(4)[1]; double reg72=0.56758792732771999991*e.pos(2)[0]; double reg73=0.56758792732771999991*e.pos(2)[1];
   double reg74=1.78379396366386*e.pos(3)[1]; double reg75=reg2+reg56; double reg76=pow(reg50,2); reg17=reg17-reg57; double reg77=pow(reg31,2);
   reg48=reg48-reg6; reg47=reg47-reg45; double reg78=1.78379396366386*e.pos(3)[0]; double reg79=reg5+reg59; double reg80=reg1+reg73;
   double reg81=0.56758792732771999991*e.pos(2)[2]; reg76=reg67+reg76; reg67=pow(reg3,2); reg63=reg68+reg63; reg68=reg4+reg72;
   reg64=reg71+reg64; double reg82=3.2673902919218320001*e.pos(5)[1]; reg47=reg61+reg47; reg66=reg39-reg66; reg6=reg36-reg6;
   reg36=0.63369514596091600003*e.pos(2)[1]; reg39=0.63369514596091600003*e.pos(2)[0]; double reg83=pow(reg48,2); double reg84=1.78379396366386*e.pos(3)[2]; double reg85=reg13+reg24;
   double reg86=3.2673902919218320001*e.pos(4)[2]; reg75=reg75-reg74; reg17=reg70+reg17; double reg87=3.2673902919218320001*e.pos(5)[0]; double reg88=0.43241207267228000009*e.pos(3)[0];
   reg4=reg59-reg4; reg14=reg14-reg58; reg62=reg62-reg69; reg1=reg56-reg1; reg56=pow(reg54,2);
   reg60=reg77+reg60; reg59=0.43241207267228000009*e.pos(3)[1]; reg79=reg79-reg78; reg77=3.2673902919218320001*e.pos(5)[2]; double reg89=0.63369514596091600003*e.pos(2)[2];
   reg69=reg47-reg69; reg47=0.366304854039084*e.pos(3)[1]; double reg90=pow(reg66,2); reg22=reg22+reg36; reg14=reg86+reg14;
   double reg91=0.366304854039084*e.pos(3)[0]; reg20=reg20+reg39; reg82=reg64-reg82; reg83=reg56+reg83; reg56=pow(reg6,2);
   reg76=reg67+reg76; reg64=pow(reg62,2); reg78=reg78+reg68; reg74=reg74+reg80; reg87=reg17-reg87;
   reg17=reg10+reg81; reg59=reg1-reg59; reg75=reg19+reg75; reg1=1.3513818909915799999*e.pos(5)[1]; reg10=reg24-reg10;
   reg88=reg4-reg88; reg85=reg85-reg84; reg4=1.3513818909915799999*e.pos(5)[0]; reg24=0.43241207267228000009*e.pos(3)[2]; reg79=reg25+reg79;
   reg60=pow(reg60,0.5); reg63=pow(reg63,0.5); reg83=reg64+reg83; reg64=reg37/reg63; reg84=reg84+reg17;
   reg67=reg43/reg63; reg77=reg14-reg77; reg19=reg19-reg74; reg79=reg79-reg4; reg14=pow(reg87,2);
   double reg92=reg27/reg60; double reg93=reg20+reg91; reg25=reg25-reg78; reg88=reg29+reg88; reg29=reg28/reg60;
   double reg94=pow(reg82,2); double reg95=3.2673902919218320001*e.pos(3)[1]; double reg96=reg35-reg36; double reg97=2.2673902919218320001*e.pos(2)[0]; reg24=reg10-reg24;
   reg10=3.2673902919218320001*e.pos(3)[0]; double reg98=reg34-reg39; reg75=reg75-reg1; reg56=reg90+reg56; reg85=reg38+reg85;
   reg90=1.3513818909915799999*e.pos(5)[2]; reg76=pow(reg76,0.5); double reg99=2.2673902919218320001*e.pos(2)[1]; double reg100=pow(reg69,2); double reg101=0.366304854039084*e.pos(3)[2];
   reg30=reg30+reg89; reg59=reg41+reg59; reg41=reg22+reg47; double reg102=reg29*reg59; reg60=reg31/reg60;
   reg55=reg24+reg55; reg24=reg92*reg88; reg94=reg14+reg94; reg14=3.2673902919218320001*e.pos(3)[2]; double reg103=reg46-reg89;
   reg95=reg96-reg95; reg96=reg67*reg75; reg10=reg98-reg10; reg56=reg100+reg56; reg85=reg85-reg90;
   reg98=reg30+reg101; reg100=reg52-reg41; double reg104=2.9010854378827480001*e.pos(5)[1]; double reg105=2.9010854378827480001*e.pos(5)[0]; double reg106=reg53-reg93;
   double reg107=reg64*reg79; reg83=pow(reg83,0.5); reg63=reg9/reg63; double reg108=pow(reg77,2); double reg109=2.2673902919218320001*e.pos(2)[2];
   reg19=reg1+reg19; reg97=reg34+reg97; reg25=reg4+reg25; reg38=reg38-reg84; reg99=reg35+reg99;
   double reg110=reg49/reg76; double reg111=reg50/reg76; double reg112=5.42101086242752217e-20*e.pos(5)[1]; double reg113=5.42101086242752217e-20*e.pos(5)[0]; reg10=reg70+reg10;
   reg76=reg3/reg76; reg70=reg63*reg85; reg96=reg107+reg96; reg14=reg103-reg14; reg103=reg110*reg25;
   reg38=reg90+reg38; reg56=pow(reg56,0.5); reg107=reg111*reg19; double reg114=reg54/reg83; double reg115=reg48/reg83;
   reg47=reg99-reg47; reg91=reg97-reg91; reg94=reg108+reg94; reg97=2.9010854378827480001*e.pos(5)[2]; reg99=reg61-reg98;
   reg106=reg106+reg105; reg95=reg71+reg95; reg100=reg104+reg100; reg102=reg24+reg102; reg24=reg60*reg55;
   reg109=reg46+reg109; reg91=reg53+reg91; reg94=pow(reg94,0.5); reg107=reg103+reg107; reg47=reg52+reg47;
   reg52=reg66/reg56; reg53=5.42101086242752217e-20*e.pos(5)[2]; reg71=reg115*reg100; reg103=reg6/reg56; reg86=reg14+reg86;
   reg99=reg99+reg97; reg14=reg114*reg106; reg108=reg76*reg38; reg83=reg62/reg83; reg10=reg10-reg113;
   reg101=reg109-reg101; reg95=reg95-reg112; reg70=reg96+reg70; reg24=reg102+reg24; reg96=reg64*reg70;
   reg86=reg86-reg53; reg102=reg67*reg70; reg108=reg107+reg108; reg71=reg14+reg71; reg47=reg47-reg104;
   reg56=reg69/reg56; reg91=reg91-reg105; reg14=reg52*reg10; reg101=reg61+reg101; reg61=reg92*reg24;
   reg107=reg103*reg95; reg109=reg29*reg24; double reg116=reg83*reg99; double reg117=reg87/reg94; double reg118=reg82/reg94;
   double reg119=reg63*reg70; double reg120=reg117*reg91; reg116=reg71+reg116; reg102=reg75-reg102; reg96=reg79-reg96;
   reg71=reg60*reg24; double reg121=reg110*reg108; double reg122=reg111*reg108; double reg123=reg118*reg47; double reg124=reg56*reg86;
   reg61=reg88-reg61; reg107=reg14+reg107; reg109=reg59-reg109; reg101=reg101-reg97; reg94=reg77/reg94;
   reg119=reg85-reg119; reg124=reg107+reg124; reg14=pow(reg102,2); reg107=pow(reg96,2); reg123=reg120+reg123;
   reg120=reg94*reg101; double reg125=reg114*reg116; reg71=reg55-reg71; reg121=reg25-reg121; double reg126=pow(reg109,2);
   reg122=reg19-reg122; double reg127=pow(reg61,2); double reg128=reg115*reg116; double reg129=reg76*reg108; reg14=reg107+reg14;
   reg107=pow(reg119,2); double reg130=pow(reg71,2); double reg131=pow(reg121,2); reg125=reg106-reg125; double reg132=reg52*reg124;
   double reg133=pow(reg122,2); double reg134=reg103*reg124; reg128=reg100-reg128; reg129=reg38-reg129; double reg135=reg83*reg116;
   reg120=reg123+reg120; reg127=reg126+reg127; reg123=reg118*reg120; reg130=reg127+reg130; reg126=pow(reg128,2);
   reg107=reg14+reg107; reg135=reg99-reg135; reg14=reg117*reg120; reg127=pow(reg125,2); double reg136=pow(reg129,2);
   reg132=reg10-reg132; double reg137=reg56*reg124; reg133=reg131+reg133; reg134=reg95-reg134; reg136=reg133+reg136;
   reg14=reg91-reg14; reg130=pow(reg130,0.5); reg126=reg127+reg126; reg127=pow(reg135,2); reg137=reg86-reg137;
   reg131=pow(reg134,2); reg133=pow(reg132,2); double reg138=reg94*reg120; reg107=pow(reg107,0.5); reg123=reg47-reg123;
   reg109=reg109/reg130; reg127=reg126+reg127; reg126=pow(reg123,2); reg136=pow(reg136,0.5); double reg139=pow(reg137,2);
   reg131=reg133+reg131; reg61=reg61/reg130; reg138=reg101-reg138; reg102=reg102/reg107; reg133=pow(reg14,2);
   reg96=reg96/reg107; reg59=reg59*reg109; reg79=reg79*reg96; reg127=pow(reg127,0.5); reg130=reg71/reg130;
   reg122=reg122/reg136; reg121=reg121/reg136; reg139=reg131+reg139; reg88=reg88*reg61; reg96=reg37*reg96;
   reg71=reg43*reg102; reg107=reg119/reg107; reg119=pow(reg138,2); reg64=reg37*reg64; reg109=reg28*reg109;
   reg61=reg27*reg61; reg126=reg133+reg126; reg29=reg28*reg29; reg102=reg75*reg102; reg92=reg27*reg92;
   reg67=reg43*reg67; reg63=reg9*reg63; reg139=pow(reg139,0.5); reg67=reg64+reg67; reg109=reg61+reg109;
   reg125=reg125/reg127; reg27=reg31*reg130; reg60=reg31*reg60; reg29=reg92+reg29; reg102=reg79+reg102;
   reg130=reg55*reg130; reg59=reg88+reg59; reg128=reg128/reg127; reg85=reg85*reg107; reg25=reg25*reg121;
   reg19=reg19*reg122; reg110=reg49*reg110; reg111=reg50*reg111; reg122=reg50*reg122; reg119=reg126+reg119;
   reg121=reg49*reg121; reg136=reg129/reg136; reg71=reg96+reg71; reg107=reg9*reg107; reg85=reg102+reg85;
   reg111=reg110+reg111; reg115=reg48*reg115; reg114=reg54*reg114; reg107=reg71+reg107; reg63=reg67+reg63;
   reg130=reg59+reg130; reg76=reg3*reg76; reg38=reg38*reg136; reg119=pow(reg119,0.5); reg122=reg121+reg122;
   reg136=reg3*reg136; reg48=reg48*reg128; reg54=reg54*reg125; reg19=reg25+reg19; reg127=reg135/reg127;
   reg128=reg100*reg128; reg125=reg106*reg125; reg132=reg132/reg139; reg134=reg134/reg139; reg60=reg29+reg60;
   reg27=reg109+reg27; reg136=reg122+reg136; reg107=reg70*reg107; reg123=reg123/reg119; reg76=reg111+reg76;
   reg38=reg19+reg38; reg103=reg6*reg103; reg52=reg66*reg52; reg3=reg62*reg127; reg48=reg54+reg48;
   reg14=reg14/reg119; reg27=reg24*reg27; reg127=reg99*reg127; reg128=reg125+reg128; reg10=reg10*reg132;
   reg95=reg95*reg134; reg139=reg137/reg139; reg132=reg66*reg132; reg134=reg6*reg134; reg130=reg60*reg130;
   reg83=reg62*reg83; reg115=reg114+reg115; reg85=reg63*reg85; reg27=reg130-reg27; reg119=reg138/reg119;
   reg56=reg69*reg56; reg103=reg52+reg103; reg3=reg48+reg3; reg127=reg128+reg127; reg95=reg10+reg95;
   reg86=reg86*reg139; reg6=reg82*reg123; reg134=reg132+reg134; reg139=reg69*reg139; reg9=reg87*reg14;
   reg83=reg115+reg83; reg107=reg85-reg107; reg14=reg91*reg14; reg123=reg47*reg123; reg136=reg108*reg136;
   reg117=reg87*reg117; reg38=reg76*reg38; reg118=reg82*reg118; reg123=reg14+reg123; reg101=reg101*reg119;
   reg10=0.005384432036113586778*reg27; reg136=reg38-reg136; reg14=0.088847818743090689935*reg27; reg19=0.005384432036113586778*reg107; reg24=0.009463616120767210603*reg107;
   reg127=reg83*reg127; reg3=reg116*reg3; reg25=0.009463616120767210603*reg27; reg119=reg77*reg119; reg28=0.088847818743090689935*reg107;
   reg6=reg9+reg6; reg56=reg103+reg56; reg86=reg95+reg86; reg139=reg134+reg139; reg9=0.021537728144454347112*reg27;
   reg29=0.021537728144454347112*reg107; reg118=reg117+reg118; reg94=reg77*reg94; reg25=reg19+reg25; reg19=reg10+reg19;
   reg31=0.009463616120767210603*reg136; reg37=reg29+reg9; reg38=0.088847818743090689935*reg136; reg28=reg9+reg28; reg9=0.021537728144454347112*reg136;
   reg14=reg29+reg14; reg94=reg118+reg94; reg86=reg56*reg86; reg139=reg124*reg139; reg119=reg6+reg119;
   reg3=reg127-reg3; reg101=reg123+reg101; reg10=reg24+reg10; reg6=0.005384432036113586778*reg136; reg10=reg10+reg6;
   reg28=reg28+reg9; reg24=0.028457289286966203713*reg3; reg139=reg86-reg139; reg38=reg37+reg38; reg29=0.016449618187943419918*reg3;
   reg37=0.0018441552587796664112*reg3; reg31=reg19+reg31; reg19=0.0041124045469858549794*reg3; reg14=reg9+reg14; reg101=reg94*reg101;
   reg25=reg6+reg25; reg119=reg120*reg119; reg37=reg28+reg37; reg6=0.016449618187943419918*reg139; reg14=reg29+reg14;
   reg9=0.016449618187943419916*reg139; reg28=0.0018441552587796664109*reg139; reg43=0.028457289286966203713*reg139; reg25=reg25+reg19; reg119=reg101-reg119;
   reg31=reg19+reg31; reg19=0.0041124045469858549794*reg139; reg47=0.004112404546985854979*reg139; reg10=reg24-reg10; reg29=reg38+reg29;
   reg47=reg10-reg47; reg28=reg14+reg28; reg19=reg31+reg19; reg10=0.028457289286966203713*reg119; reg14=0.0018441552587796664111*reg119;
   reg6=reg37+reg6; reg24=0.016449618187943419918*reg119; reg31=0.0041124045469858549794*reg119; reg25=reg43-reg25; reg9=reg29+reg9;
   Ne(0,15)+=reg24+reg28; Ne(1,16)+=reg24+reg28; Ne(2,17)+=reg24+reg28; Ne(0,0)+=reg47-reg31; Ne(1,1)+=reg47-reg31;
   Ne(2,2)+=reg47-reg31; Ne(0,12)+=reg6+reg24; Ne(1,13)+=reg6+reg24; Ne(2,14)+=reg6+reg24; Ne(0,3)+=reg25-reg31;
   Ne(1,4)+=reg25-reg31; Ne(2,5)+=reg25-reg31; Ne(0,6)+=reg10-reg19; Ne(1,7)+=reg10-reg19; Ne(2,8)+=reg10-reg19;
   Ne(0,9)+=reg9+reg14; Ne(1,10)+=reg9+reg14; Ne(2,11)+=reg9+reg14;

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
   double reg70=reg60-reg53; reg38=reg38+reg48; double reg71=1.8213672050459180516*e.pos(4)[1]; reg18=reg18-reg60; double reg72=0.66666666666666670528*e.pos(5)[1];
   double reg73=0.12200846792814621817*e.pos(3)[0]; double reg74=0.36602540378443865451*e.pos(3)[2]; double reg75=reg54-reg55; reg49=reg50+reg49; double reg76=reg68+reg67;
   reg9=reg9+reg69; double reg77=1.8213672050459180516*e.pos(6)[0]; reg75=reg69+reg75; reg18=reg18+reg72; double reg78=reg62-reg65;
   double reg79=1.8213672050459180516*e.pos(6)[1]; reg70=reg72+reg70; reg38=reg38-reg71; reg37=reg37-reg64; double reg80=reg73+reg63;
   reg49=reg49+reg74; double reg81=reg66+reg61; double reg82=1.8213672050459180516*e.pos(4)[2]; double reg83=0.66666666666666670528*e.pos(5)[2]; reg27=reg27-reg62;
   double reg84=0.12200846792814621817*e.pos(3)[2]; double reg85=0.48803387171258487271*e.pos(6)[1]; reg38=reg72+reg38; reg18=reg18-reg79; double reg86=0.66666666666666670528*e.pos(7)[1];
   reg75=reg77+reg75; double reg87=reg64-reg80; reg70=reg79+reg70; double reg88=reg3+reg25; reg9=reg9-reg77;
   double reg89=0.66666666666666670528*e.pos(7)[0]; double reg90=reg4+reg26; double reg91=0.48803387171258487271*e.pos(6)[0]; reg37=reg69+reg37; reg78=reg83+reg78;
   double reg92=1.8213672050459180516*e.pos(6)[2]; double reg93=reg0+reg21; double reg94=reg84+reg76; reg27=reg27+reg83; double reg95=reg7+reg20;
   reg49=reg49-reg82; double reg96=reg71-reg81; reg19=reg5+reg19; reg5=0.48803387171258487271*e.pos(6)[2]; reg38=reg38-reg85;
   reg75=reg75-reg89; reg49=reg83+reg49; reg17=reg6+reg17; reg87=reg69+reg87; reg70=reg70-reg86;
   reg37=reg37-reg91; reg78=reg92+reg78; reg34=reg34+reg93; reg35=reg35+reg95; reg6=reg23+reg36;
   reg69=reg82-reg94; reg96=reg72+reg96; reg72=reg15+reg45; reg16=reg16+reg90; reg9=reg9-reg89;
   double reg97=0.66666666666666670528*e.pos(7)[2]; reg27=reg27-reg92; reg8=reg8+reg88; reg18=reg18-reg86; reg38=reg38-reg86;
   reg41=reg1+reg41; reg73=reg8+reg73; reg37=reg37-reg89; reg1=0.66666666666666670528*e.pos(4)[0]; reg49=reg49-reg5;
   reg34=reg29+reg34; reg66=reg16+reg66; reg8=0.66666666666666670528*e.pos(4)[1]; reg78=reg78-reg97; reg35=reg31+reg35;
   reg87=reg91+reg87; reg19=reg13+reg19; reg16=pow(reg18,2); reg29=pow(reg9,2); reg31=pow(reg75,2);
   reg40=reg2+reg40; reg69=reg83+reg69; reg27=reg27-reg97; reg50=reg50+reg6; reg39=reg22+reg39;
   reg2=pow(reg70,2); reg96=reg85+reg96; reg17=reg10+reg17; reg24=reg24+reg72; reg86=reg96-reg86;
   reg84=reg24+reg84; reg22=0.66666666666666670528*e.pos(4)[2]; reg66=reg66-reg8; reg24=1.8213672050459180516*e.pos(5)[1]; reg50=reg51+reg50;
   reg51=pow(reg37,2); reg83=0.48803387171258487271*e.pos(5)[0]; reg35=reg35-reg8; reg96=0.48803387171258487271*e.pos(5)[1]; reg69=reg5+reg69;
   reg34=reg34-reg1; reg44=reg14+reg44; reg2=reg31+reg2; reg14=pow(reg27,2); reg31=reg48+reg17;
   double reg98=reg47+reg19; reg39=reg33+reg39; reg89=reg87-reg89; reg87=pow(reg78,2); reg49=reg49-reg97;
   reg40=reg46+reg40; reg41=reg56+reg41; reg16=reg29+reg16; reg29=pow(reg38,2); reg73=reg73-reg1;
   double reg99=1.8213672050459180516*e.pos(5)[0]; double reg100=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg22; double reg101=reg74+reg39; double reg102=reg8+reg31;
   reg97=reg69-reg97; reg69=pow(reg86,2); double reg103=reg1+reg98; double reg104=0.66666666666666670528*e.pos(6)[1]; reg66=reg66-reg24;
   reg14=reg16+reg14; reg16=pow(reg89,2); double reg105=pow(reg49,2); reg29=reg51+reg29; reg51=0.66666666666666670528*e.pos(6)[0];
   reg73=reg73-reg99; reg50=reg50-reg22; double reg106=0.48803387171258487271*e.pos(5)[2]; reg44=reg67+reg44; reg35=reg35-reg96;
   reg34=reg34-reg83; double reg107=reg42+reg41; reg87=reg2+reg87; reg2=reg43+reg40; reg87=pow(reg87,0.5);
   reg69=reg16+reg69; reg16=1.8213672050459180516*e.pos(7)[1]; double reg108=reg99-reg103; double reg109=reg58+reg44; double reg110=0.48803387171258487271*e.pos(7)[1];
   reg66=reg66+reg104; reg35=reg104+reg35; reg14=pow(reg14,0.5); double reg111=1.8213672050459180516*e.pos(7)[0]; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg73=reg73+reg51; reg34=reg51+reg34; reg105=reg29+reg105; reg29=reg24-reg102; double reg113=reg22+reg101;
   double reg114=pow(reg97,2); reg8=reg8+reg2; reg50=reg50-reg106; double reg115=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg100;
   reg1=reg1+reg107; double reg116=reg9/reg14; double reg117=reg75/reg87; double reg118=reg18/reg14; double reg119=reg83-reg1;
   double reg120=reg100-reg113; double reg121=reg70/reg87; double reg122=reg96-reg8; reg73=reg73-reg112; reg34=reg34-reg111;
   reg105=pow(reg105,0.5); reg84=reg84+reg115; double reg123=0.48803387171258487271*e.pos(7)[2]; reg108=reg51+reg108; reg22=reg22+reg109;
   reg29=reg104+reg29; reg114=reg69+reg114; reg66=reg66-reg110; reg35=reg35-reg16; reg69=1.8213672050459180516*e.pos(7)[2];
   reg50=reg115+reg50; reg14=reg27/reg14; reg108=reg112+reg108; reg119=reg51+reg119; reg120=reg115+reg120;
   reg122=reg104+reg122; reg51=reg106-reg22; reg104=reg38/reg105; double reg124=reg37/reg105; reg29=reg110+reg29;
   reg50=reg50-reg69; reg114=pow(reg114,0.5); reg84=reg84-reg123; double reg125=reg121*reg35; double reg126=reg118*reg66;
   double reg127=reg117*reg34; double reg128=reg116*reg73; reg87=reg78/reg87; double reg129=reg104*reg29; reg120=reg123+reg120;
   double reg130=reg124*reg108; reg125=reg127+reg125; reg126=reg128+reg126; reg122=reg16+reg122; reg127=reg14*reg84;
   reg51=reg115+reg51; reg115=reg89/reg114; reg105=reg49/reg105; reg128=reg86/reg114; double reg131=reg87*reg50;
   reg119=reg111+reg119; double reg132=reg105*reg120; double reg133=reg128*reg122; reg131=reg125+reg131; reg129=reg130+reg129;
   reg114=reg97/reg114; reg125=reg115*reg119; reg127=reg126+reg127; reg51=reg69+reg51; reg126=reg118*reg127;
   reg130=reg116*reg127; reg132=reg129+reg132; reg133=reg125+reg133; reg125=reg114*reg51; reg129=reg121*reg131;
   double reg134=reg117*reg131; reg129=reg35-reg129; double reg135=reg124*reg132; reg134=reg34-reg134; double reg136=reg104*reg132;
   reg125=reg133+reg125; reg133=reg87*reg131; reg126=reg66-reg126; reg130=reg73-reg130; double reg137=reg14*reg127;
   double reg138=pow(reg126,2); double reg139=pow(reg130,2); double reg140=reg105*reg132; reg136=reg29-reg136; reg137=reg84-reg137;
   double reg141=reg128*reg125; reg135=reg108-reg135; double reg142=reg115*reg125; reg133=reg50-reg133; double reg143=pow(reg134,2);
   double reg144=pow(reg129,2); reg142=reg119-reg142; double reg145=reg114*reg125; reg141=reg122-reg141; double reg146=pow(reg135,2);
   double reg147=pow(reg136,2); reg140=reg120-reg140; double reg148=pow(reg137,2); reg138=reg139+reg138; reg139=pow(reg133,2);
   reg144=reg143+reg144; reg147=reg146+reg147; reg143=pow(reg140,2); reg139=reg144+reg139; reg148=reg138+reg148;
   reg145=reg51-reg145; reg138=pow(reg142,2); reg144=pow(reg141,2); reg148=pow(reg148,0.5); reg144=reg138+reg144;
   reg139=pow(reg139,0.5); reg143=reg147+reg143; reg138=pow(reg145,2); reg130=reg130/reg148; reg126=reg126/reg148;
   reg134=reg134/reg139; reg143=pow(reg143,0.5); reg129=reg129/reg139; reg138=reg144+reg138; reg34=reg34*reg134;
   reg35=reg35*reg129; reg139=reg133/reg139; reg121=reg70*reg121; reg117=reg75*reg117; reg133=reg18*reg126;
   reg144=reg9*reg130; reg148=reg137/reg148; reg126=reg66*reg126; reg130=reg73*reg130; reg118=reg18*reg118;
   reg116=reg9*reg116; reg138=pow(reg138,0.5); reg135=reg135/reg143; reg136=reg136/reg143; reg129=reg70*reg129;
   reg134=reg75*reg134; reg84=reg84*reg148; reg126=reg130+reg126; reg133=reg144+reg133; reg148=reg27*reg148;
   reg9=reg38*reg136; reg18=reg37*reg135; reg14=reg27*reg14; reg118=reg116+reg118; reg104=reg38*reg104;
   reg121=reg117+reg121; reg143=reg140/reg143; reg87=reg78*reg87; reg136=reg29*reg136; reg50=reg50*reg139;
   reg35=reg34+reg35; reg129=reg134+reg129; reg139=reg78*reg139; reg141=reg141/reg138; reg124=reg37*reg124;
   reg142=reg142/reg138; reg135=reg108*reg135; reg119=reg119*reg142; reg120=reg120*reg143; reg122=reg122*reg141;
   reg138=reg145/reg138; reg14=reg118+reg14; reg142=reg89*reg142; reg136=reg135+reg136; reg9=reg18+reg9;
   reg143=reg49*reg143; reg141=reg86*reg141; reg50=reg35+reg50; reg139=reg129+reg139; reg115=reg89*reg115;
   reg87=reg121+reg87; reg128=reg86*reg128; reg104=reg124+reg104; reg148=reg133+reg148; reg105=reg49*reg105;
   reg84=reg126+reg84; reg122=reg119+reg122; reg139=reg131*reg139; reg84=reg14*reg84; reg51=reg51*reg138;
   reg50=reg87*reg50; reg128=reg115+reg128; reg114=reg97*reg114; reg141=reg142+reg141; reg138=reg97*reg138;
   reg120=reg136+reg120; reg148=reg127*reg148; reg143=reg9+reg143; reg105=reg104+reg105; reg138=reg141+reg138;
   reg51=reg122+reg51; reg139=reg50-reg139; reg148=reg84-reg148; reg114=reg128+reg114; reg143=reg132*reg143;
   reg120=reg105*reg120; reg9=0.024056261216234395431*reg148; reg14=0.13144585576580215187*reg139; reg18=0.024056261216234395431*reg139; reg27=0.024056261216234409915*reg139;
   reg29=0.13144585576580215187*reg148; reg34=0.04166666666666666908*reg148; reg35=0.035220810900864524453*reg148; reg37=0.024056261216234409915*reg148; reg38=0.035220810900864524453*reg139;
   reg138=reg125*reg138; reg51=reg114*reg51; reg143=reg120-reg143; reg49=0.04166666666666666908*reg139; reg50=reg35+reg38;
   reg66=0.024056261216234409915*reg143; reg18=reg18-reg34; reg70=0.13144585576580215187*reg143; reg38=reg38+reg29; reg9=reg9-reg49;
   reg29=reg29+reg14; reg73=0.035220810900864524453*reg143; reg27=reg34+reg27; reg34=0.024056261216234395431*reg143; reg14=reg35+reg14;
   reg138=reg51-reg138; reg35=0.04166666666666666908*reg143; reg49=reg37+reg49; reg37=0.035220810900864524453*reg138; reg38=reg70+reg38;
   reg51=0.13144585576580215187*reg138; reg29=reg29+reg73; reg14=reg73+reg14; reg70=reg50+reg70; reg66=reg18-reg66;
   reg18=0.024056261216234409915*reg138; reg9=reg9-reg35; reg50=0.04166666666666666908*reg138; reg27=reg34-reg27; reg34=0.024056261216234395431*reg138;
   reg35=reg49+reg35; Ne(0,9)+=reg66-reg50; Ne(1,10)+=reg66-reg50; Ne(2,11)+=reg66-reg50; Ne(0,12)+=reg70+reg51;
   Ne(1,13)+=reg70+reg51; Ne(2,14)+=reg70+reg51; Ne(0,6)+=reg9-reg18; Ne(1,7)+=reg9-reg18; Ne(2,8)+=reg9-reg18;
   Ne(0,15)+=reg38+reg37; Ne(1,16)+=reg38+reg37; Ne(2,17)+=reg38+reg37; Ne(0,3)+=reg27-reg50; Ne(1,4)+=reg27-reg50;
   Ne(2,5)+=reg27-reg50; Ne(0,0)+=reg34-reg35; Ne(1,1)+=reg34-reg35; Ne(2,2)+=reg34-reg35; Ne(0,18)+=reg37+reg29;
   Ne(1,19)+=reg37+reg29; Ne(2,20)+=reg37+reg29; Ne(0,21)+=reg51+reg14; Ne(1,22)+=reg51+reg14; Ne(2,23)+=reg51+reg14;

}

};
} // namespace LMT

