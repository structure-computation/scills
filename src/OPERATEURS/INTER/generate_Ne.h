
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
   
   double reg0=0.56758792732771999991*e.pos(1)[1]; double reg1=0.78379396366385999995*e.pos(0)[1]; double reg2=0.78379396366385999995*e.pos(1)[1]; double reg3=0.56758792732771999991*e.pos(0)[1]; double reg4=0.78379396366385999995*e.pos(1)[0];
   double reg5=0.56758792732771999991*e.pos(0)[0]; double reg6=0.78379396366385999995*e.pos(0)[0]; double reg7=0.56758792732771999991*e.pos(1)[0]; double reg8=reg4+reg5; double reg9=1.3513818909915799999*e.pos(3)[0];
   double reg10=reg0+reg1; double reg11=0.78379396366385999995*e.pos(0)[2]; double reg12=0.56758792732771999991*e.pos(1)[2]; double reg13=0.78379396366385999995*e.pos(1)[2]; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg7+reg6; double reg16=reg2+reg3; double reg17=1.3513818909915799999*e.pos(3)[1]; reg16=reg16-reg17; double reg18=1.78379396366386*e.pos(4)[1];
   double reg19=reg17-reg10; double reg20=2.2673902919218320001*e.pos(0)[1]; double reg21=reg13+reg14; double reg22=1.3513818909915799999*e.pos(3)[2]; double reg23=reg9-reg15;
   double reg24=2.2673902919218320001*e.pos(0)[0]; double reg25=0.63369514596091600003*e.pos(1)[1]; double reg26=reg12+reg11; double reg27=1.78379396366386*e.pos(4)[0]; reg8=reg8-reg9;
   double reg28=0.63369514596091600003*e.pos(1)[0]; double reg29=0.63369514596091600003*e.pos(0)[1]; double reg30=2.2673902919218320001*e.pos(1)[1]; double reg31=0.63369514596091600003*e.pos(0)[0]; double reg32=2.2673902919218320001*e.pos(1)[0];
   double reg33=2.2673902919218320001*e.pos(0)[2]; double reg34=0.63369514596091600003*e.pos(1)[2]; double reg35=2.9010854378827480001*e.pos(3)[1]; double reg36=reg25+reg20; double reg37=reg28+reg24;
   reg23=reg27+reg23; reg19=reg18+reg19; double reg38=reg22-reg26; double reg39=2.9010854378827480001*e.pos(3)[0]; double reg40=0.43241207267228000009*e.pos(4)[0];
   double reg41=0.43241207267228000009*e.pos(4)[1]; reg2=reg2-reg1; reg4=reg4-reg6; double reg42=1.78379396366386*e.pos(5)[0]; reg8=reg8+reg27;
   double reg43=1.78379396366386*e.pos(4)[2]; reg21=reg21-reg22; reg16=reg16+reg18; double reg44=1.78379396366386*e.pos(5)[1]; reg32=reg32+reg31;
   reg8=reg8-reg42; reg19=reg19-reg44; reg44=reg16-reg44; reg16=2.9010854378827480001*e.pos(3)[2]; double reg45=reg34+reg33;
   reg42=reg23-reg42; reg21=reg21+reg43; reg23=1.78379396366386*e.pos(5)[2]; double reg46=0.366304854039084*e.pos(4)[1]; double reg47=reg35-reg36;
   reg4=reg40+reg4; double reg48=reg39-reg37; double reg49=0.366304854039084*e.pos(4)[0]; double reg50=0.43241207267228000009*e.pos(5)[0]; double reg51=0.43241207267228000009*e.pos(4)[2];
   reg2=reg41+reg2; double reg52=0.43241207267228000009*e.pos(5)[1]; double reg53=0.63369514596091600003*e.pos(0)[2]; reg13=reg13-reg11; double reg54=2.2673902919218320001*e.pos(1)[2];
   reg30=reg30+reg29; reg38=reg43+reg38; double reg55=pow(reg44,2); reg13=reg51+reg13; double reg56=pow(reg42,2);
   double reg57=0.43241207267228000009*e.pos(5)[2]; reg52=reg2-reg52; reg2=pow(reg8,2); double reg58=pow(reg19,2); reg50=reg4-reg50;
   reg21=reg21-reg23; reg23=reg38-reg23; reg48=reg48+reg49; reg4=0.366304854039084*e.pos(5)[0]; reg32=reg32-reg39;
   reg38=2.710505431213761085e-20*e.pos(3)[0]; reg30=reg30-reg35; double reg59=0.366304854039084*e.pos(4)[2]; double reg60=reg16-reg45; double reg61=reg31-reg28;
   double reg62=0.78379396366385999995*e.pos(2)[1]; reg54=reg54+reg53; double reg63=0.366304854039084*e.pos(5)[1]; reg47=reg47+reg46; double reg64=reg29-reg25;
   double reg65=0.78379396366385999995*e.pos(2)[0]; double reg66=2.710505431213761085e-20*e.pos(3)[1]; double reg67=reg53-reg34; double reg68=0.78379396366385999995*e.pos(2)[2]; double reg69=2.710505431213761085e-20*e.pos(3)[2];
   double reg70=3.2673902919218320001*e.pos(4)[1]; double reg71=0.56758792732771999991*e.pos(2)[1]; double reg72=reg5+reg65; double reg73=0.56758792732771999991*e.pos(2)[0]; reg57=reg13-reg57;
   reg30=reg46+reg30; reg13=1.78379396366386*e.pos(3)[1]; double reg74=pow(reg52,2); double reg75=1.78379396366386*e.pos(3)[0]; double reg76=3.2673902919218320001*e.pos(4)[0];
   double reg77=reg3+reg62; double reg78=pow(reg50,2); reg54=reg54-reg16; double reg79=0.366304854039084*e.pos(5)[2]; reg60=reg60+reg59;
   reg61=reg61-reg38; reg55=reg2+reg55; reg47=reg47-reg63; reg58=reg56+reg58; reg32=reg49+reg32;
   reg64=reg64-reg66; reg2=pow(reg21,2); reg56=pow(reg23,2); reg48=reg48-reg4; reg61=reg76+reg61;
   double reg80=3.2673902919218320001*e.pos(5)[0]; reg54=reg59+reg54; double reg81=pow(reg47,2); double reg82=pow(reg48,2); reg67=reg67-reg69;
   reg74=reg78+reg74; reg63=reg30-reg63; reg60=reg60-reg79; reg30=0.63369514596091600003*e.pos(2)[0]; reg78=0.63369514596091600003*e.pos(2)[1];
   double reg83=pow(reg57,2); double reg84=reg6+reg73; double reg85=reg1+reg71; double reg86=0.56758792732771999991*e.pos(2)[2]; reg4=reg32-reg4;
   reg32=3.2673902919218320001*e.pos(4)[2]; double reg87=3.2673902919218320001*e.pos(5)[1]; reg64=reg70+reg64; reg56=reg58+reg56; reg72=reg72-reg75;
   reg1=reg62-reg1; reg58=0.43241207267228000009*e.pos(3)[1]; reg6=reg65-reg6; reg77=reg77-reg13; reg62=0.43241207267228000009*e.pos(3)[0];
   reg65=1.78379396366386*e.pos(3)[2]; reg2=reg55+reg2; reg55=reg14+reg68; reg13=reg13+reg85; reg55=reg55-reg65;
   reg87=reg64-reg87; reg64=pow(reg4,2); reg83=reg74+reg83; reg74=0.43241207267228000009*e.pos(3)[2]; reg68=reg68-reg11;
   reg62=reg6-reg62; reg75=reg75+reg84; reg58=reg1-reg58; reg1=1.3513818909915799999*e.pos(5)[1]; reg6=pow(reg60,2);
   reg80=reg61-reg80; reg61=0.63369514596091600003*e.pos(2)[2]; reg20=reg20+reg78; double reg88=pow(reg63,2); reg79=reg54-reg79;
   reg2=pow(reg2,0.5); reg11=reg11+reg86; reg81=reg82+reg81; reg72=reg27+reg72; reg54=0.366304854039084*e.pos(3)[0];
   reg24=reg24+reg30; reg82=0.366304854039084*e.pos(3)[1]; double reg89=1.3513818909915799999*e.pos(5)[0]; reg77=reg18+reg77; reg56=pow(reg56,0.5);
   double reg90=3.2673902919218320001*e.pos(5)[2]; reg67=reg32+reg67; double reg91=reg44/reg2; reg27=reg27-reg75; double reg92=reg42/reg56;
   reg88=reg64+reg88; reg64=reg19/reg56; reg18=reg18-reg13; reg72=reg72-reg89; reg65=reg65+reg11;
   double reg93=pow(reg79,2); double reg94=reg31-reg30; double reg95=3.2673902919218320001*e.pos(3)[0]; double reg96=reg29-reg78; double reg97=3.2673902919218320001*e.pos(3)[1];
   double reg98=reg20+reg82; reg6=reg81+reg6; reg33=reg33+reg61; reg81=2.2673902919218320001*e.pos(2)[1]; double reg99=reg24+reg54;
   double reg100=0.366304854039084*e.pos(3)[2]; double reg101=2.2673902919218320001*e.pos(2)[0]; double reg102=1.3513818909915799999*e.pos(5)[2]; reg90=reg67-reg90; reg67=reg8/reg2;
   double reg103=pow(reg80,2); reg83=pow(reg83,0.5); reg40=reg62+reg40; reg41=reg58+reg41; reg77=reg77-reg1;
   reg74=reg68-reg74; reg58=pow(reg87,2); reg55=reg43+reg55; reg62=reg33+reg100; reg68=reg91*reg77;
   reg55=reg55-reg102; reg93=reg88+reg93; reg88=reg67*reg72; double reg104=2.9010854378827480001*e.pos(5)[1]; double reg105=reg46-reg98;
   reg95=reg94-reg95; reg6=pow(reg6,0.5); reg43=reg43-reg65; reg18=reg1+reg18; reg56=reg23/reg56;
   reg27=reg89+reg27; reg2=reg21/reg2; reg94=reg52/reg83; double reg106=reg50/reg83; double reg107=reg92*reg40;
   double reg108=reg64*reg41; reg51=reg74+reg51; reg58=reg103+reg58; reg74=pow(reg90,2); reg101=reg31+reg101;
   reg81=reg29+reg81; reg103=reg49-reg99; double reg109=2.2673902919218320001*e.pos(2)[2]; double reg110=2.9010854378827480001*e.pos(5)[0]; reg97=reg96-reg97;
   reg96=reg53-reg61; double reg111=3.2673902919218320001*e.pos(3)[2]; double reg112=reg56*reg51; reg103=reg103+reg110; double reg113=5.42101086242752217e-20*e.pos(5)[0];
   reg74=reg58+reg74; reg58=reg2*reg55; double reg114=2.9010854378827480001*e.pos(5)[2]; double reg115=reg59-reg62; reg54=reg101-reg54;
   reg82=reg81-reg82; reg81=reg47/reg6; reg105=reg105+reg104; reg101=5.42101086242752217e-20*e.pos(5)[1]; reg70=reg97+reg70;
   reg109=reg53+reg109; reg97=reg48/reg6; reg43=reg102+reg43; double reg116=reg94*reg18; reg93=pow(reg93,0.5);
   double reg117=reg106*reg27; reg83=reg57/reg83; reg111=reg96-reg111; reg76=reg95+reg76; reg108=reg107+reg108;
   reg68=reg88+reg68; reg88=reg83*reg43; reg100=reg109-reg100; reg95=5.42101086242752217e-20*e.pos(5)[2]; reg116=reg117+reg116;
   reg82=reg46+reg82; reg112=reg108+reg112; reg46=reg81*reg105; reg96=reg4/reg93; reg107=reg97*reg103;
   reg6=reg60/reg6; reg54=reg49+reg54; reg70=reg70-reg101; reg76=reg76-reg113; reg49=reg63/reg93;
   reg32=reg111+reg32; reg115=reg115+reg114; reg74=pow(reg74,0.5); reg58=reg68+reg58; reg46=reg107+reg46;
   reg68=reg6*reg115; reg93=reg79/reg93; reg107=reg96*reg76; reg108=reg49*reg70; reg32=reg32-reg95;
   reg109=reg64*reg112; reg111=reg80/reg74; reg117=reg87/reg74; reg54=reg54-reg110; reg82=reg82-reg104;
   reg100=reg59+reg100; reg59=reg67*reg58; double reg118=reg92*reg112; double reg119=reg91*reg58; reg88=reg116+reg88;
   reg116=reg117*reg82; double reg120=reg56*reg112; reg68=reg46+reg68; reg46=reg111*reg54; reg74=reg90/reg74;
   reg109=reg41-reg109; reg118=reg40-reg118; double reg121=reg2*reg58; double reg122=reg93*reg32; reg119=reg77-reg119;
   double reg123=reg106*reg88; reg59=reg72-reg59; reg108=reg107+reg108; reg107=reg94*reg88; reg100=reg100-reg114;
   double reg124=pow(reg118,2); reg116=reg46+reg116; reg46=reg74*reg100; reg120=reg51-reg120; double reg125=pow(reg109,2);
   reg122=reg108+reg122; reg123=reg27-reg123; reg107=reg18-reg107; reg108=reg83*reg88; reg121=reg55-reg121;
   double reg126=pow(reg59,2); double reg127=pow(reg119,2); double reg128=reg81*reg68; double reg129=reg97*reg68; double reg130=pow(reg121,2);
   double reg131=pow(reg123,2); double reg132=pow(reg107,2); double reg133=pow(reg120,2); reg108=reg43-reg108; reg129=reg103-reg129;
   reg127=reg126+reg127; reg124=reg125+reg124; reg128=reg105-reg128; reg125=reg49*reg122; reg126=reg6*reg68;
   reg46=reg116+reg46; reg116=reg96*reg122; reg125=reg70-reg125; double reg134=reg93*reg122; double reg135=reg117*reg46;
   reg116=reg76-reg116; reg130=reg127+reg130; reg132=reg131+reg132; reg133=reg124+reg133; reg124=pow(reg129,2);
   reg127=pow(reg128,2); reg131=pow(reg108,2); double reg136=reg111*reg46; reg126=reg115-reg126; reg135=reg82-reg135;
   double reg137=reg74*reg46; reg131=reg132+reg131; reg132=pow(reg116,2); double reg138=pow(reg125,2); reg134=reg32-reg134;
   reg130=pow(reg130,0.5); reg133=pow(reg133,0.5); reg136=reg54-reg136; reg127=reg124+reg127; reg124=pow(reg126,2);
   reg138=reg132+reg138; reg132=pow(reg134,2); reg109=reg109/reg133; double reg139=pow(reg135,2); reg118=reg118/reg133;
   double reg140=pow(reg136,2); reg137=reg100-reg137; reg131=pow(reg131,0.5); reg124=reg127+reg124; reg59=reg59/reg130;
   reg119=reg119/reg130; reg133=reg120/reg133; reg64=reg19*reg64; reg120=reg42*reg118; reg92=reg42*reg92;
   reg19=reg19*reg109; reg42=reg44*reg119; reg127=reg8*reg59; reg119=reg77*reg119; reg132=reg138+reg132;
   reg77=pow(reg137,2); reg139=reg140+reg139; reg130=reg121/reg130; reg59=reg72*reg59; reg124=pow(reg124,0.5);
   reg91=reg44*reg91; reg67=reg8*reg67; reg107=reg107/reg131; reg118=reg40*reg118; reg123=reg123/reg131;
   reg109=reg41*reg109; reg131=reg108/reg131; reg128=reg128/reg124; reg109=reg118+reg109; reg129=reg129/reg124;
   reg132=pow(reg132,0.5); reg18=reg18*reg107; reg27=reg27*reg123; reg119=reg59+reg119; reg91=reg67+reg91;
   reg8=reg23*reg133; reg42=reg127+reg42; reg40=reg130*reg21; reg55=reg130*reg55; reg19=reg120+reg19;
   reg77=reg139+reg77; reg56=reg23*reg56; reg64=reg92+reg64; reg2=reg21*reg2; reg123=reg50*reg123;
   reg107=reg52*reg107; reg133=reg51*reg133; reg94=reg52*reg94; reg106=reg50*reg106; reg103=reg103*reg129;
   reg43=reg43*reg131; reg81=reg47*reg81; reg97=reg48*reg97; reg47=reg47*reg128; reg56=reg64+reg56;
   reg94=reg106+reg94; reg83=reg57*reg83; reg2=reg91+reg2; reg133=reg109+reg133; reg131=reg57*reg131;
   reg107=reg123+reg107; reg77=pow(reg77,0.5); reg40=reg42+reg40; reg129=reg48*reg129; reg119=reg55+reg119;
   reg116=reg116/reg132; reg8=reg19+reg8; reg125=reg125/reg132; reg18=reg27+reg18; reg124=reg126/reg124;
   reg128=reg105*reg128; reg19=reg60*reg124; reg47=reg129+reg47; reg8=reg112*reg8; reg124=reg115*reg124;
   reg133=reg56*reg133; reg128=reg103+reg128; reg49=reg63*reg49; reg43=reg18+reg43; reg96=reg4*reg96;
   reg131=reg107+reg131; reg135=reg135/reg77; reg136=reg136/reg77; reg40=reg58*reg40; reg119=reg2*reg119;
   reg83=reg94+reg83; reg81=reg97+reg81; reg63=reg63*reg125; reg4=reg4*reg116; reg6=reg60*reg6;
   reg132=reg134/reg132; reg125=reg70*reg125; reg116=reg76*reg116; reg54=reg54*reg136; reg82=reg82*reg135;
   reg77=reg137/reg77; reg124=reg128+reg124; reg6=reg81+reg6; reg136=reg80*reg136; reg135=reg87*reg135;
   reg19=reg47+reg19; reg93=reg79*reg93; reg131=reg88*reg131; reg40=reg119-reg40; reg125=reg116+reg125;
   reg32=reg32*reg132; reg63=reg4+reg63; reg132=reg79*reg132; reg49=reg96+reg49; reg8=reg133-reg8;
   reg43=reg83*reg43; reg111=reg80*reg111; reg117=reg87*reg117; reg2=reg90*reg77; reg135=reg136+reg135;
   reg93=reg49+reg93; reg4=0.088847818743090689935*reg40; reg77=reg100*reg77; reg19=reg68*reg19; reg82=reg54+reg82;
   reg131=reg43-reg131; reg18=0.005384432036113586778*reg40; reg21=0.009463616120767210603*reg40; reg74=reg90*reg74; reg23=0.021537728144454347112*reg40;
   reg117=reg111+reg117; reg32=reg125+reg32; reg124=reg6*reg124; reg132=reg63+reg132; reg6=0.009463616120767210603*reg8;
   reg27=0.088847818743090689935*reg8; reg41=0.021537728144454347112*reg8; reg42=0.005384432036113586778*reg8; reg43=0.009463616120767210603*reg131; reg6=reg18+reg6;
   reg44=0.088847818743090689935*reg131; reg4=reg41+reg4; reg41=reg23+reg41; reg18=reg42+reg18; reg19=reg124-reg19;
   reg74=reg117+reg74; reg42=reg21+reg42; reg132=reg122*reg132; reg32=reg93*reg32; reg27=reg23+reg27;
   reg21=0.005384432036113586778*reg131; reg77=reg82+reg77; reg23=0.021537728144454347112*reg131; reg2=reg135+reg2; reg47=0.016449618187943419918*reg19;
   reg48=0.0018441552587796664112*reg19; reg44=reg41+reg44; reg4=reg4+reg23; reg27=reg23+reg27; reg42=reg42+reg21;
   reg132=reg32-reg132; reg43=reg18+reg43; reg18=0.0041124045469858549794*reg19; reg77=reg74*reg77; reg6=reg21+reg6;
   reg2=reg46*reg2; reg21=0.028457289286966203713*reg19; reg42=reg21-reg42; reg27=reg47+reg27; reg21=0.016449618187943419918*reg132;
   reg48=reg4+reg48; reg4=0.004112404546985854979*reg132; reg2=reg77-reg2; reg23=0.0041124045469858549794*reg132; reg32=0.016449618187943419916*reg132;
   reg47=reg44+reg47; reg41=0.0018441552587796664109*reg132; reg44=0.028457289286966203713*reg132; reg6=reg6+reg18; reg43=reg18+reg43;
   reg23=reg43+reg23; reg18=0.028457289286966203713*reg2; reg43=0.016449618187943419918*reg2; reg21=reg48+reg21; reg41=reg27+reg41;
   reg6=reg44-reg6; reg4=reg42-reg4; reg27=0.0041124045469858549794*reg2; reg42=0.0018441552587796664111*reg2; reg32=reg47+reg32;
   Ne(0,15)+=reg43+reg41; Ne(1,16)+=reg43+reg41; Ne(2,17)+=reg43+reg41; Ne(0,12)+=reg21+reg43; Ne(1,13)+=reg21+reg43;
   Ne(2,14)+=reg21+reg43; Ne(0,9)+=reg32+reg42; Ne(1,10)+=reg32+reg42; Ne(2,11)+=reg32+reg42; Ne(0,0)+=reg4-reg27;
   Ne(1,1)+=reg4-reg27; Ne(2,2)+=reg4-reg27; Ne(0,3)+=reg6-reg27; Ne(1,4)+=reg6-reg27; Ne(2,5)+=reg6-reg27;
   Ne(0,6)+=reg18-reg23; Ne(1,7)+=reg18-reg23; Ne(2,8)+=reg18-reg23;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(1)[0]; double reg1=0.12200846792814621817*e.pos(1)[0]; double reg2=0.12200846792814621817*e.pos(1)[1]; double reg3=0.36602540378443865451*e.pos(0)[0]; double reg4=0.36602540378443865451*e.pos(0)[1];
   double reg5=0.12200846792814621817*e.pos(0)[0]; double reg6=0.12200846792814621817*e.pos(0)[1]; double reg7=0.36602540378443865451*e.pos(1)[1]; double reg8=reg1+reg3; double reg9=1.3660254037844385386*e.pos(2)[0];
   double reg10=0.45534180126147951289*e.pos(2)[1]; double reg11=reg7+reg6; double reg12=reg0+reg5; double reg13=0.45534180126147951289*e.pos(2)[0]; double reg14=0.12200846792814621817*e.pos(1)[2];
   double reg15=0.36602540378443865451*e.pos(0)[2]; double reg16=1.3660254037844385386*e.pos(2)[1]; double reg17=1.3660254037844385386*e.pos(1)[1]; double reg18=reg2+reg4; double reg19=1.3660254037844385386*e.pos(1)[0];
   double reg20=0.45534180126147951289*e.pos(0)[1]; double reg21=0.45534180126147951289*e.pos(0)[0]; double reg22=0.12200846792814621817*e.pos(0)[2]; double reg23=0.36602540378443865451*e.pos(1)[2]; double reg24=reg14+reg15;
   double reg25=1.3660254037844385386*e.pos(2)[2]; double reg26=0.45534180126147951289*e.pos(1)[0]; double reg27=0.45534180126147951289*e.pos(1)[1]; double reg28=reg12+reg13; double reg29=1.3660254037844385386*e.pos(3)[0];
   double reg30=reg11+reg10; double reg31=1.3660254037844385386*e.pos(3)[1]; double reg32=reg23+reg22; double reg33=0.45534180126147951289*e.pos(2)[2]; double reg34=0.12200846792814621817*e.pos(2)[0];
   double reg35=0.12200846792814621817*e.pos(2)[1]; double reg36=0.45534180126147951289*e.pos(0)[2]; double reg37=reg21+reg19; double reg38=reg20+reg17; double reg39=1.3660254037844385386*e.pos(1)[2];
   double reg40=1.3660254037844385386*e.pos(0)[1]; double reg41=1.3660254037844385386*e.pos(0)[0]; double reg42=0.45534180126147951289*e.pos(3)[0]; reg8=reg8+reg9; double reg43=0.45534180126147951289*e.pos(3)[1];
   reg18=reg18+reg16; double reg44=reg36+reg39; double reg45=0.36602540378443865451*e.pos(3)[1]; reg38=reg35+reg38; double reg46=0.36602540378443865451*e.pos(3)[0];
   reg37=reg34+reg37; double reg47=0.36602540378443865451*e.pos(2)[1]; double reg48=0.45534180126147951289*e.pos(1)[2]; double reg49=1.3660254037844385386*e.pos(0)[2]; double reg50=0.12200846792814621817*e.pos(2)[2];
   double reg51=1.3660254037844385386*e.pos(3)[2]; double reg52=reg32+reg33; double reg53=reg30+reg31; reg8=reg8+reg42; double reg54=0.48803387171258487271*e.pos(4)[0];
   double reg55=reg28+reg29; double reg56=0.45534180126147951289*e.pos(3)[2]; reg24=reg24+reg25; reg18=reg18+reg43; double reg57=0.36602540378443865451*e.pos(2)[0];
   double reg58=reg26+reg41; double reg59=reg27+reg40; double reg60=0.48803387171258487271*e.pos(4)[1]; double reg61=reg48+reg49; double reg62=0.36602540378443865451*e.pos(2)[2];
   reg24=reg24+reg56; reg8=reg8-reg54; double reg63=reg59+reg47; double reg64=0.48803387171258487271*e.pos(4)[2]; double reg65=reg58+reg57;
   reg37=reg37+reg46; double reg66=1.8213672050459180516*e.pos(4)[0]; double reg67=0.12200846792814621817*e.pos(3)[1]; double reg68=reg52+reg51; double reg69=0.66666666666666670528*e.pos(5)[0];
   reg44=reg50+reg44; double reg70=0.36602540378443865451*e.pos(3)[2]; double reg71=reg54-reg55; double reg72=reg60-reg53; reg18=reg18-reg60;
   reg38=reg38+reg45; double reg73=1.8213672050459180516*e.pos(4)[1]; double reg74=0.12200846792814621817*e.pos(3)[0]; double reg75=0.66666666666666670528*e.pos(5)[1]; double reg76=1.8213672050459180516*e.pos(6)[0];
   double reg77=reg61+reg62; reg71=reg69+reg71; reg18=reg18+reg75; double reg78=1.8213672050459180516*e.pos(6)[1]; double reg79=reg64-reg68;
   reg72=reg75+reg72; reg8=reg8+reg69; double reg80=1.8213672050459180516*e.pos(4)[2]; reg44=reg44+reg70; reg38=reg38-reg73;
   double reg81=reg74+reg65; reg37=reg37-reg66; double reg82=0.66666666666666670528*e.pos(5)[2]; reg24=reg24-reg64; double reg83=reg67+reg63;
   double reg84=0.12200846792814621817*e.pos(3)[2]; reg44=reg44-reg80; reg8=reg8-reg76; double reg85=reg3+reg26; reg71=reg76+reg71;
   double reg86=0.48803387171258487271*e.pos(6)[1]; reg38=reg75+reg38; double reg87=reg73-reg83; reg18=reg18-reg78; double reg88=0.66666666666666670528*e.pos(7)[1];
   double reg89=reg66-reg81; reg72=reg78+reg72; double reg90=reg84+reg77; double reg91=0.66666666666666670528*e.pos(7)[0]; double reg92=reg4+reg27;
   double reg93=0.48803387171258487271*e.pos(6)[0]; reg37=reg69+reg37; double reg94=reg7+reg20; reg79=reg82+reg79; double reg95=1.8213672050459180516*e.pos(6)[2];
   double reg96=reg0+reg21; reg24=reg24+reg82; double reg97=0.48803387171258487271*e.pos(6)[2]; reg44=reg82+reg44; reg35=reg35+reg94;
   reg87=reg75+reg87; reg19=reg5+reg19; reg5=reg23+reg36; reg38=reg38-reg86; reg34=reg34+reg96;
   reg71=reg71-reg91; reg75=reg80-reg90; reg79=reg95+reg79; reg17=reg6+reg17; reg89=reg69+reg89;
   reg72=reg72-reg88; reg37=reg37-reg93; reg6=reg15+reg48; reg16=reg16+reg92; reg8=reg8-reg91;
   reg69=0.66666666666666670528*e.pos(7)[2]; reg24=reg24-reg95; reg9=reg9+reg85; reg18=reg18-reg88; reg38=reg38-reg88;
   reg41=reg1+reg41; reg74=reg9+reg74; reg37=reg37-reg91; reg1=0.66666666666666670528*e.pos(4)[0]; reg44=reg44-reg97;
   reg34=reg29+reg34; reg67=reg16+reg67; reg9=0.66666666666666670528*e.pos(4)[1]; reg79=reg79-reg69; reg35=reg31+reg35;
   reg89=reg93+reg89; reg19=reg13+reg19; reg16=pow(reg18,2); reg29=pow(reg8,2); reg31=pow(reg71,2);
   reg40=reg2+reg40; reg75=reg82+reg75; reg24=reg24-reg69; reg50=reg50+reg5; reg39=reg22+reg39;
   reg2=pow(reg72,2); reg87=reg86+reg87; reg17=reg10+reg17; reg25=reg25+reg6; reg88=reg87-reg88;
   reg84=reg25+reg84; reg22=0.66666666666666670528*e.pos(4)[2]; reg67=reg67-reg9; reg25=1.8213672050459180516*e.pos(5)[1]; reg50=reg51+reg50;
   reg51=pow(reg37,2); reg82=0.48803387171258487271*e.pos(5)[0]; reg35=reg35-reg9; reg87=0.48803387171258487271*e.pos(5)[1]; reg75=reg97+reg75;
   reg34=reg34-reg1; reg49=reg14+reg49; reg2=reg31+reg2; reg14=pow(reg24,2); reg31=reg45+reg17;
   double reg98=reg46+reg19; reg39=reg33+reg39; reg91=reg89-reg91; reg89=pow(reg79,2); reg44=reg44-reg69;
   reg40=reg47+reg40; reg41=reg57+reg41; reg16=reg29+reg16; reg29=pow(reg38,2); reg74=reg74-reg1;
   double reg99=1.8213672050459180516*e.pos(5)[0]; double reg100=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg22; double reg101=reg70+reg39; double reg102=reg9+reg31;
   reg69=reg75-reg69; reg75=pow(reg88,2); double reg103=reg1+reg98; double reg104=0.66666666666666670528*e.pos(6)[1]; reg67=reg67-reg25;
   reg14=reg16+reg14; reg16=pow(reg91,2); double reg105=pow(reg44,2); reg29=reg51+reg29; reg51=0.66666666666666670528*e.pos(6)[0];
   reg74=reg74-reg99; reg50=reg50-reg22; double reg106=0.48803387171258487271*e.pos(5)[2]; reg49=reg62+reg49; reg35=reg35-reg87;
   reg34=reg34-reg82; double reg107=reg42+reg41; reg89=reg2+reg89; reg2=reg43+reg40; reg89=pow(reg89,0.5);
   reg75=reg16+reg75; reg16=1.8213672050459180516*e.pos(7)[1]; double reg108=reg99-reg103; double reg109=reg56+reg49; double reg110=0.48803387171258487271*e.pos(7)[1];
   reg67=reg67+reg104; reg35=reg104+reg35; reg14=pow(reg14,0.5); double reg111=1.8213672050459180516*e.pos(7)[0]; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg74=reg74+reg51; reg34=reg51+reg34; reg105=reg29+reg105; reg29=reg25-reg102; double reg113=reg22+reg101;
   double reg114=pow(reg69,2); reg9=reg9+reg2; reg50=reg50-reg106; double reg115=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg100;
   reg1=reg1+reg107; double reg116=reg8/reg14; double reg117=reg71/reg89; double reg118=reg18/reg14; double reg119=reg82-reg1;
   double reg120=reg100-reg113; double reg121=reg72/reg89; double reg122=reg87-reg9; reg74=reg74-reg112; reg34=reg34-reg111;
   reg105=pow(reg105,0.5); reg84=reg84+reg115; double reg123=0.48803387171258487271*e.pos(7)[2]; reg108=reg51+reg108; reg22=reg22+reg109;
   reg29=reg104+reg29; reg114=reg75+reg114; reg67=reg67-reg110; reg35=reg35-reg16; reg75=1.8213672050459180516*e.pos(7)[2];
   reg50=reg115+reg50; reg14=reg24/reg14; reg108=reg112+reg108; reg119=reg51+reg119; reg120=reg115+reg120;
   reg122=reg104+reg122; reg51=reg106-reg22; reg104=reg38/reg105; double reg124=reg37/reg105; reg29=reg110+reg29;
   reg50=reg50-reg75; reg114=pow(reg114,0.5); reg84=reg84-reg123; double reg125=reg121*reg35; double reg126=reg118*reg67;
   double reg127=reg117*reg34; double reg128=reg116*reg74; reg89=reg79/reg89; double reg129=reg104*reg29; reg120=reg123+reg120;
   double reg130=reg124*reg108; reg125=reg127+reg125; reg126=reg128+reg126; reg122=reg16+reg122; reg127=reg14*reg84;
   reg51=reg115+reg51; reg115=reg91/reg114; reg105=reg44/reg105; reg128=reg88/reg114; double reg131=reg89*reg50;
   reg119=reg111+reg119; double reg132=reg105*reg120; double reg133=reg128*reg122; reg131=reg125+reg131; reg129=reg130+reg129;
   reg114=reg69/reg114; reg125=reg115*reg119; reg127=reg126+reg127; reg51=reg75+reg51; reg126=reg118*reg127;
   reg130=reg116*reg127; reg132=reg129+reg132; reg133=reg125+reg133; reg125=reg114*reg51; reg129=reg121*reg131;
   double reg134=reg117*reg131; reg129=reg35-reg129; double reg135=reg124*reg132; reg134=reg34-reg134; double reg136=reg104*reg132;
   reg125=reg133+reg125; reg133=reg89*reg131; reg126=reg67-reg126; reg130=reg74-reg130; double reg137=reg14*reg127;
   double reg138=pow(reg126,2); double reg139=pow(reg130,2); double reg140=reg105*reg132; reg136=reg29-reg136; reg137=reg84-reg137;
   double reg141=reg128*reg125; reg135=reg108-reg135; double reg142=reg115*reg125; reg133=reg50-reg133; double reg143=pow(reg134,2);
   double reg144=pow(reg129,2); reg142=reg119-reg142; double reg145=reg114*reg125; reg141=reg122-reg141; double reg146=pow(reg135,2);
   double reg147=pow(reg136,2); reg140=reg120-reg140; double reg148=pow(reg137,2); reg138=reg139+reg138; reg139=pow(reg133,2);
   reg144=reg143+reg144; reg147=reg146+reg147; reg143=pow(reg140,2); reg139=reg144+reg139; reg148=reg138+reg148;
   reg145=reg51-reg145; reg138=pow(reg142,2); reg144=pow(reg141,2); reg148=pow(reg148,0.5); reg144=reg138+reg144;
   reg139=pow(reg139,0.5); reg143=reg147+reg143; reg138=pow(reg145,2); reg130=reg130/reg148; reg126=reg126/reg148;
   reg134=reg134/reg139; reg143=pow(reg143,0.5); reg129=reg129/reg139; reg138=reg144+reg138; reg34=reg34*reg134;
   reg35=reg35*reg129; reg139=reg133/reg139; reg121=reg72*reg121; reg117=reg71*reg117; reg133=reg18*reg126;
   reg144=reg8*reg130; reg148=reg137/reg148; reg126=reg67*reg126; reg130=reg74*reg130; reg118=reg18*reg118;
   reg116=reg8*reg116; reg138=pow(reg138,0.5); reg135=reg135/reg143; reg136=reg136/reg143; reg129=reg72*reg129;
   reg134=reg71*reg134; reg84=reg84*reg148; reg126=reg130+reg126; reg133=reg144+reg133; reg148=reg24*reg148;
   reg8=reg38*reg136; reg18=reg37*reg135; reg14=reg24*reg14; reg118=reg116+reg118; reg104=reg38*reg104;
   reg121=reg117+reg121; reg143=reg140/reg143; reg89=reg79*reg89; reg136=reg29*reg136; reg50=reg50*reg139;
   reg35=reg34+reg35; reg129=reg134+reg129; reg139=reg79*reg139; reg141=reg141/reg138; reg124=reg37*reg124;
   reg142=reg142/reg138; reg135=reg108*reg135; reg119=reg119*reg142; reg120=reg120*reg143; reg122=reg122*reg141;
   reg138=reg145/reg138; reg14=reg118+reg14; reg142=reg91*reg142; reg136=reg135+reg136; reg8=reg18+reg8;
   reg143=reg44*reg143; reg141=reg88*reg141; reg50=reg35+reg50; reg139=reg129+reg139; reg115=reg91*reg115;
   reg89=reg121+reg89; reg128=reg88*reg128; reg104=reg124+reg104; reg148=reg133+reg148; reg105=reg44*reg105;
   reg84=reg126+reg84; reg122=reg119+reg122; reg139=reg131*reg139; reg84=reg14*reg84; reg51=reg51*reg138;
   reg50=reg89*reg50; reg128=reg115+reg128; reg114=reg69*reg114; reg141=reg142+reg141; reg138=reg69*reg138;
   reg120=reg136+reg120; reg148=reg127*reg148; reg143=reg8+reg143; reg105=reg104+reg105; reg138=reg141+reg138;
   reg51=reg122+reg51; reg139=reg50-reg139; reg148=reg84-reg148; reg114=reg128+reg114; reg143=reg132*reg143;
   reg120=reg105*reg120; reg8=0.024056261216234395431*reg148; reg14=0.13144585576580215187*reg139; reg18=0.024056261216234395431*reg139; reg24=0.024056261216234409915*reg139;
   reg29=0.13144585576580215187*reg148; reg34=0.04166666666666666908*reg148; reg35=0.035220810900864524453*reg148; reg37=0.024056261216234409915*reg148; reg38=0.035220810900864524453*reg139;
   reg138=reg125*reg138; reg51=reg114*reg51; reg143=reg120-reg143; reg44=0.04166666666666666908*reg139; reg50=reg35+reg38;
   reg67=0.024056261216234409915*reg143; reg18=reg18-reg34; reg69=0.13144585576580215187*reg143; reg38=reg38+reg29; reg8=reg8-reg44;
   reg29=reg29+reg14; reg71=0.035220810900864524453*reg143; reg24=reg34+reg24; reg34=0.024056261216234395431*reg143; reg14=reg35+reg14;
   reg138=reg51-reg138; reg35=0.04166666666666666908*reg143; reg44=reg37+reg44; reg37=0.035220810900864524453*reg138; reg38=reg69+reg38;
   reg51=0.13144585576580215187*reg138; reg29=reg29+reg71; reg14=reg71+reg14; reg69=reg50+reg69; reg67=reg18-reg67;
   reg18=0.024056261216234409915*reg138; reg8=reg8-reg35; reg50=0.04166666666666666908*reg138; reg24=reg34-reg24; reg34=0.024056261216234395431*reg138;
   reg35=reg44+reg35; Ne(0,9)+=reg67-reg50; Ne(1,10)+=reg67-reg50; Ne(2,11)+=reg67-reg50; Ne(0,12)+=reg69+reg51;
   Ne(1,13)+=reg69+reg51; Ne(2,14)+=reg69+reg51; Ne(0,6)+=reg8-reg18; Ne(1,7)+=reg8-reg18; Ne(2,8)+=reg8-reg18;
   Ne(0,15)+=reg38+reg37; Ne(1,16)+=reg38+reg37; Ne(2,17)+=reg38+reg37; Ne(0,3)+=reg24-reg50; Ne(1,4)+=reg24-reg50;
   Ne(2,5)+=reg24-reg50; Ne(0,0)+=reg34-reg35; Ne(1,1)+=reg34-reg35; Ne(2,2)+=reg34-reg35; Ne(0,18)+=reg37+reg29;
   Ne(1,19)+=reg37+reg29; Ne(2,20)+=reg37+reg29; Ne(0,21)+=reg51+reg14; Ne(1,22)+=reg51+reg14; Ne(2,23)+=reg51+reg14;

}

};
} // namespace LMT

