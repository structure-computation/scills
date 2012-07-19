
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
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.56758792732771999991*e.pos(0)[1]; double reg3=0.78379396366385999995*e.pos(0)[1]; double reg4=0.56758792732771999991*e.pos(1)[0];
   double reg5=0.56758792732771999991*e.pos(1)[1]; double reg6=0.78379396366385999995*e.pos(0)[0]; double reg7=0.78379396366385999995*e.pos(1)[1]; double reg8=reg7+reg2; double reg9=1.3513818909915799999*e.pos(3)[1];
   double reg10=1.3513818909915799999*e.pos(3)[0]; double reg11=reg0+reg1; double reg12=reg5+reg3; double reg13=0.78379396366385999995*e.pos(1)[2]; double reg14=0.56758792732771999991*e.pos(1)[2];
   double reg15=0.78379396366385999995*e.pos(0)[2]; double reg16=0.56758792732771999991*e.pos(0)[2]; double reg17=reg4+reg6; double reg18=0.63369514596091600003*e.pos(1)[1]; reg8=reg8-reg9;
   double reg19=2.2673902919218320001*e.pos(0)[1]; double reg20=reg10-reg17; double reg21=1.78379396366386*e.pos(4)[1]; double reg22=1.3513818909915799999*e.pos(3)[2]; double reg23=reg13+reg16;
   double reg24=reg9-reg12; double reg25=1.78379396366386*e.pos(4)[0]; reg11=reg11-reg10; double reg26=0.63369514596091600003*e.pos(1)[0]; double reg27=2.2673902919218320001*e.pos(0)[0];
   double reg28=reg14+reg15; double reg29=2.2673902919218320001*e.pos(0)[2]; reg20=reg25+reg20; double reg30=reg18+reg19; double reg31=0.63369514596091600003*e.pos(1)[2];
   reg7=reg7-reg3; reg0=reg0-reg6; double reg32=0.43241207267228000009*e.pos(4)[1]; double reg33=2.9010854378827480001*e.pos(3)[1]; double reg34=0.63369514596091600003*e.pos(0)[1];
   double reg35=2.2673902919218320001*e.pos(1)[1]; double reg36=0.43241207267228000009*e.pos(4)[0]; double reg37=0.63369514596091600003*e.pos(0)[0]; double reg38=2.2673902919218320001*e.pos(1)[0]; double reg39=reg22-reg28;
   double reg40=1.78379396366386*e.pos(4)[2]; reg23=reg23-reg22; double reg41=reg26+reg27; double reg42=1.78379396366386*e.pos(5)[1]; reg8=reg8+reg21;
   reg24=reg21+reg24; double reg43=2.9010854378827480001*e.pos(3)[0]; double reg44=1.78379396366386*e.pos(5)[0]; reg11=reg11+reg25; reg7=reg32+reg7;
   double reg45=0.366304854039084*e.pos(4)[1]; double reg46=reg33-reg30; double reg47=reg43-reg41; double reg48=0.366304854039084*e.pos(4)[0]; reg13=reg13-reg15;
   double reg49=0.43241207267228000009*e.pos(5)[1]; reg39=reg40+reg39; reg24=reg24-reg42; reg11=reg11-reg44; reg42=reg8-reg42;
   reg23=reg23+reg40; reg8=1.78379396366386*e.pos(5)[2]; double reg50=reg31+reg29; double reg51=2.9010854378827480001*e.pos(3)[2]; double reg52=0.43241207267228000009*e.pos(5)[0];
   double reg53=0.43241207267228000009*e.pos(4)[2]; reg38=reg38+reg37; reg0=reg36+reg0; reg35=reg35+reg34; double reg54=2.2673902919218320001*e.pos(1)[2];
   double reg55=0.63369514596091600003*e.pos(0)[2]; reg44=reg20-reg44; reg13=reg53+reg13; reg20=pow(reg24,2); reg49=reg7-reg49;
   reg52=reg0-reg52; reg0=pow(reg44,2); reg39=reg39-reg8; reg7=0.43241207267228000009*e.pos(5)[2]; reg8=reg23-reg8;
   reg23=reg37-reg26; reg54=reg54+reg55; double reg56=pow(reg42,2); reg35=reg35-reg33; double reg57=0.78379396366385999995*e.pos(2)[0];
   double reg58=pow(reg11,2); reg38=reg38-reg43; double reg59=2.710505431213761085e-20*e.pos(3)[0]; double reg60=0.366304854039084*e.pos(4)[2]; double reg61=0.78379396366385999995*e.pos(2)[1];
   double reg62=reg51-reg50; double reg63=reg34-reg18; double reg64=0.366304854039084*e.pos(5)[1]; reg46=reg46+reg45; double reg65=2.710505431213761085e-20*e.pos(3)[1];
   double reg66=0.366304854039084*e.pos(5)[0]; reg47=reg47+reg48; double reg67=pow(reg39,2); reg23=reg23-reg59; reg63=reg63-reg65;
   double reg68=reg55-reg31; double reg69=2.710505431213761085e-20*e.pos(3)[2]; double reg70=3.2673902919218320001*e.pos(4)[1]; double reg71=3.2673902919218320001*e.pos(4)[0]; reg54=reg54-reg51;
   reg35=reg45+reg35; reg38=reg48+reg38; double reg72=0.366304854039084*e.pos(5)[2]; reg62=reg62+reg60; reg46=reg46-reg64;
   double reg73=pow(reg52,2); reg47=reg47-reg66; double reg74=pow(reg49,2); double reg75=0.56758792732771999991*e.pos(2)[1]; reg7=reg13-reg7;
   reg13=0.56758792732771999991*e.pos(2)[0]; double reg76=0.78379396366385999995*e.pos(2)[2]; double reg77=1.78379396366386*e.pos(3)[0]; double reg78=reg1+reg57; double reg79=pow(reg8,2);
   double reg80=reg2+reg61; reg56=reg58+reg56; reg20=reg0+reg20; reg0=1.78379396366386*e.pos(3)[1]; reg79=reg56+reg79;
   reg56=0.56758792732771999991*e.pos(2)[2]; reg78=reg78-reg77; reg58=pow(reg46,2); reg80=reg80-reg0; reg62=reg62-reg72;
   reg68=reg68-reg69; reg66=reg38-reg66; reg38=3.2673902919218320001*e.pos(4)[2]; double reg81=3.2673902919218320001*e.pos(5)[1]; reg63=reg70+reg63;
   double reg82=reg3+reg75; double reg83=3.2673902919218320001*e.pos(5)[0]; reg23=reg71+reg23; reg67=reg20+reg67; reg20=pow(reg47,2);
   double reg84=pow(reg7,2); reg54=reg60+reg54; reg57=reg57-reg6; double reg85=0.43241207267228000009*e.pos(3)[0]; reg64=reg35-reg64;
   reg35=1.78379396366386*e.pos(3)[2]; reg74=reg73+reg74; reg73=reg16+reg76; reg3=reg61-reg3; reg61=0.43241207267228000009*e.pos(3)[1];
   double reg86=0.63369514596091600003*e.pos(2)[0]; double reg87=0.63369514596091600003*e.pos(2)[1]; reg6=reg6+reg13; double reg88=1.3513818909915799999*e.pos(5)[0]; double reg89=pow(reg62,2);
   reg72=reg54-reg72; reg79=pow(reg79,0.5); reg54=pow(reg66,2); reg77=reg77+reg6; reg84=reg74+reg84;
   reg74=pow(reg64,2); reg27=reg27+reg86; double reg90=0.366304854039084*e.pos(3)[0]; reg67=pow(reg67,0.5); double reg91=0.63369514596091600003*e.pos(2)[2];
   double reg92=0.366304854039084*e.pos(3)[1]; reg73=reg73-reg35; reg19=reg19+reg87; reg85=reg57-reg85; reg78=reg25+reg78;
   reg61=reg3-reg61; reg3=reg15+reg56; reg15=reg76-reg15; reg80=reg21+reg80; reg83=reg23-reg83;
   reg0=reg0+reg82; reg81=reg63-reg81; reg23=1.3513818909915799999*e.pos(5)[1]; reg68=reg38+reg68; reg57=3.2673902919218320001*e.pos(5)[2];
   reg58=reg20+reg58; reg20=0.43241207267228000009*e.pos(3)[2]; reg25=reg25-reg77; reg63=reg27+reg90; reg21=reg21-reg0;
   reg35=reg35+reg3; reg76=reg19+reg92; double reg93=0.366304854039084*e.pos(3)[2]; reg29=reg29+reg91; reg74=reg54+reg74;
   reg73=reg40+reg73; reg54=2.2673902919218320001*e.pos(2)[1]; reg36=reg85+reg36; reg32=reg61+reg32; reg61=1.3513818909915799999*e.pos(5)[2];
   reg85=2.2673902919218320001*e.pos(2)[0]; reg20=reg15-reg20; reg80=reg80-reg23; reg57=reg68-reg57; reg15=reg24/reg67;
   reg68=pow(reg81,2); double reg94=pow(reg83,2); reg78=reg78-reg88; double reg95=3.2673902919218320001*e.pos(3)[1]; double reg96=reg44/reg67;
   double reg97=reg34-reg87; reg89=reg58+reg89; reg58=3.2673902919218320001*e.pos(3)[0]; double reg98=reg37-reg86; double reg99=pow(reg72,2);
   reg84=pow(reg84,0.5); double reg100=reg11/reg79; double reg101=reg42/reg79; reg73=reg73-reg61; double reg102=reg45-reg76;
   double reg103=reg101*reg80; double reg104=2.9010854378827480001*e.pos(5)[0]; double reg105=reg100*reg78; double reg106=reg48-reg63; reg89=pow(reg89,0.5);
   double reg107=2.2673902919218320001*e.pos(2)[2]; reg67=reg39/reg67; reg54=reg34+reg54; double reg108=reg96*reg36; reg85=reg37+reg85;
   double reg109=reg15*reg32; reg53=reg20+reg53; reg20=pow(reg57,2); reg68=reg94+reg68; reg94=3.2673902919218320001*e.pos(3)[2];
   double reg110=reg55-reg91; reg95=reg97-reg95; reg58=reg98-reg58; reg97=2.9010854378827480001*e.pos(5)[1]; reg99=reg74+reg99;
   reg74=reg29+reg93; reg98=reg52/reg84; double reg111=reg49/reg84; reg79=reg8/reg79; reg25=reg88+reg25;
   reg21=reg23+reg21; reg40=reg40-reg35; reg107=reg55+reg107; reg92=reg54-reg92; reg102=reg102+reg97;
   reg90=reg85-reg90; reg54=reg60-reg74; reg85=2.9010854378827480001*e.pos(5)[2]; reg103=reg105+reg103; reg20=reg68+reg20;
   reg106=reg106+reg104; reg94=reg110-reg94; reg68=5.42101086242752217e-20*e.pos(5)[1]; reg70=reg95+reg70; reg95=5.42101086242752217e-20*e.pos(5)[0];
   reg71=reg58+reg71; reg58=reg46/reg89; reg99=pow(reg99,0.5); reg105=reg47/reg89; reg110=reg79*reg73;
   reg109=reg108+reg109; reg84=reg7/reg84; reg108=reg98*reg25; double reg112=reg111*reg21; reg40=reg61+reg40;
   double reg113=reg67*reg53; reg93=reg107-reg93; reg89=reg62/reg89; reg92=reg45+reg92; reg90=reg48+reg90;
   reg45=reg105*reg106; reg48=reg58*reg102; reg107=reg84*reg40; reg54=reg54+reg85; reg20=pow(reg20,0.5);
   reg113=reg109+reg113; reg112=reg108+reg112; reg108=reg66/reg99; reg109=reg64/reg99; double reg114=5.42101086242752217e-20*e.pos(5)[2];
   reg38=reg94+reg38; reg70=reg70-reg68; reg71=reg71-reg95; reg110=reg103+reg110; reg94=reg108*reg71;
   reg103=reg100*reg110; double reg115=reg109*reg70; reg99=reg72/reg99; reg38=reg38-reg114; double reg116=reg101*reg110;
   double reg117=reg15*reg113; double reg118=reg96*reg113; double reg119=reg83/reg20; double reg120=reg81/reg20; double reg121=reg89*reg54;
   reg107=reg112+reg107; reg48=reg45+reg48; reg92=reg92-reg97; reg90=reg90-reg104; reg93=reg60+reg93;
   reg117=reg32-reg117; reg45=reg79*reg110; reg60=reg67*reg113; reg118=reg36-reg118; reg112=reg119*reg90;
   reg93=reg93-reg85; reg20=reg57/reg20; reg116=reg80-reg116; double reg122=reg120*reg92; double reg123=reg99*reg38;
   reg121=reg48+reg121; reg115=reg94+reg115; reg103=reg78-reg103; reg48=reg111*reg107; reg94=reg98*reg107;
   double reg124=pow(reg117,2); double reg125=pow(reg116,2); reg45=reg73-reg45; reg60=reg53-reg60; double reg126=reg84*reg107;
   double reg127=pow(reg118,2); reg48=reg21-reg48; double reg128=pow(reg103,2); reg123=reg115+reg123; reg122=reg112+reg122;
   reg112=reg20*reg93; reg94=reg25-reg94; reg115=reg105*reg121; double reg129=reg58*reg121; reg127=reg124+reg127;
   reg125=reg128+reg125; reg126=reg40-reg126; reg124=pow(reg45,2); reg128=pow(reg60,2); reg115=reg106-reg115;
   double reg130=reg109*reg123; reg112=reg122+reg112; reg122=reg89*reg121; double reg131=pow(reg94,2); double reg132=reg108*reg123;
   reg129=reg102-reg129; double reg133=pow(reg48,2); reg122=reg54-reg122; double reg134=pow(reg126,2); double reg135=reg119*reg112;
   double reg136=pow(reg129,2); double reg137=reg120*reg112; reg132=reg71-reg132; double reg138=reg99*reg123; reg130=reg70-reg130;
   double reg139=pow(reg115,2); reg127=reg128+reg127; reg124=reg125+reg124; reg133=reg131+reg133; reg138=reg38-reg138;
   reg125=pow(reg130,2); reg128=pow(reg132,2); reg134=reg133+reg134; reg131=reg20*reg112; reg135=reg90-reg135;
   reg137=reg92-reg137; reg136=reg139+reg136; reg133=pow(reg122,2); reg127=pow(reg127,0.5); reg124=pow(reg124,0.5);
   reg134=pow(reg134,0.5); reg139=pow(reg135,2); double reg140=pow(reg137,2); reg133=reg136+reg133; reg131=reg93-reg131;
   reg118=reg118/reg127; reg116=reg116/reg124; reg125=reg128+reg125; reg117=reg117/reg127; reg103=reg103/reg124;
   reg128=pow(reg138,2); reg136=pow(reg131,2); reg140=reg139+reg140; reg36=reg36*reg118; reg32=reg32*reg117;
   reg118=reg44*reg118; reg117=reg24*reg117; reg15=reg24*reg15; reg96=reg44*reg96; reg128=reg125+reg128;
   reg127=reg60/reg127; reg124=reg45/reg124; reg80=reg80*reg116; reg133=pow(reg133,0.5); reg78=reg78*reg103;
   reg94=reg94/reg134; reg100=reg11*reg100; reg101=reg42*reg101; reg116=reg42*reg116; reg48=reg48/reg134;
   reg103=reg11*reg103; reg80=reg78+reg80; reg98=reg52*reg98; reg111=reg49*reg111; reg128=pow(reg128,0.5);
   reg21=reg21*reg48; reg67=reg39*reg67; reg15=reg96+reg15; reg134=reg126/reg134; reg136=reg140+reg136;
   reg101=reg100+reg101; reg79=reg8*reg79; reg53=reg53*reg127; reg32=reg36+reg32; reg115=reg115/reg133;
   reg129=reg129/reg133; reg8=reg8*reg124; reg116=reg103+reg116; reg127=reg39*reg127; reg117=reg118+reg117;
   reg48=reg49*reg48; reg124=reg73*reg124; reg25=reg25*reg94; reg94=reg52*reg94; reg79=reg101+reg79;
   reg8=reg116+reg8; reg53=reg32+reg53; reg106=reg106*reg115; reg124=reg80+reg124; reg11=reg7*reg134;
   reg127=reg117+reg127; reg136=pow(reg136,0.5); reg67=reg15+reg67; reg21=reg25+reg21; reg111=reg98+reg111;
   reg58=reg46*reg58; reg134=reg40*reg134; reg105=reg47*reg105; reg46=reg46*reg129; reg115=reg47*reg115;
   reg129=reg102*reg129; reg130=reg130/reg128; reg133=reg122/reg133; reg132=reg132/reg128; reg48=reg94+reg48;
   reg84=reg7*reg84; reg7=reg62*reg133; reg46=reg115+reg46; reg129=reg106+reg129; reg133=reg54*reg133;
   reg11=reg48+reg11; reg127=reg113*reg127; reg124=reg79*reg124; reg109=reg64*reg109; reg108=reg66*reg108;
   reg71=reg71*reg132; reg70=reg70*reg130; reg8=reg110*reg8; reg128=reg138/reg128; reg134=reg21+reg134;
   reg132=reg66*reg132; reg130=reg64*reg130; reg58=reg105+reg58; reg89=reg62*reg89; reg137=reg137/reg136;
   reg135=reg135/reg136; reg53=reg67*reg53; reg84=reg111+reg84; reg11=reg107*reg11; reg134=reg84*reg134;
   reg89=reg58+reg89; reg127=reg53-reg127; reg99=reg72*reg99; reg109=reg108+reg109; reg70=reg71+reg70;
   reg38=reg38*reg128; reg136=reg131/reg136; reg130=reg132+reg130; reg128=reg72*reg128; reg92=reg92*reg137;
   reg90=reg90*reg135; reg120=reg81*reg120; reg119=reg83*reg119; reg8=reg124-reg8; reg133=reg129+reg133;
   reg7=reg46+reg7; reg137=reg81*reg137; reg135=reg83*reg135; reg15=0.005384432036113586778*reg8; reg92=reg90+reg92;
   reg21=0.009463616120767210603*reg127; reg24=0.088847818743090689935*reg127; reg25=reg57*reg136; reg136=reg93*reg136; reg137=reg135+reg137;
   reg32=0.088847818743090689935*reg8; reg36=0.021537728144454347112*reg127; reg39=0.021537728144454347112*reg8; reg20=reg57*reg20; reg120=reg119+reg120;
   reg40=0.009463616120767210603*reg8; reg7=reg121*reg7; reg133=reg89*reg133; reg11=reg134-reg11; reg128=reg130+reg128;
   reg42=0.005384432036113586778*reg127; reg38=reg70+reg38; reg99=reg109+reg99; reg24=reg39+reg24; reg21=reg15+reg21;
   reg15=reg42+reg15; reg44=0.009463616120767210603*reg11; reg45=0.021537728144454347112*reg11; reg32=reg36+reg32; reg46=0.005384432036113586778*reg11;
   reg36=reg39+reg36; reg42=reg40+reg42; reg39=0.088847818743090689935*reg11; reg20=reg120+reg20; reg25=reg137+reg25;
   reg7=reg133-reg7; reg128=reg123*reg128; reg38=reg99*reg38; reg136=reg92+reg136; reg40=0.016449618187943419918*reg7;
   reg39=reg36+reg39; reg128=reg38-reg128; reg42=reg42+reg46; reg136=reg20*reg136; reg44=reg15+reg44;
   reg32=reg32+reg45; reg15=0.028457289286966203713*reg7; reg20=0.0018441552587796664112*reg7; reg36=0.0041124045469858549794*reg7; reg21=reg46+reg21;
   reg25=reg112*reg25; reg24=reg45+reg24; reg24=reg40+reg24; reg20=reg32+reg20; reg32=0.0018441552587796664109*reg128;
   reg38=0.016449618187943419918*reg128; reg45=0.016449618187943419916*reg128; reg46=0.004112404546985854979*reg128; reg42=reg15-reg42; reg25=reg136-reg25;
   reg15=0.028457289286966203713*reg128; reg21=reg21+reg36; reg44=reg36+reg44; reg36=0.0041124045469858549794*reg128; reg40=reg39+reg40;
   reg39=0.028457289286966203713*reg25; reg36=reg44+reg36; reg21=reg15-reg21; reg32=reg24+reg32; reg15=0.0018441552587796664111*reg25;
   reg38=reg20+reg38; reg20=0.016449618187943419918*reg25; reg46=reg42-reg46; reg24=0.0041124045469858549794*reg25; reg45=reg40+reg45;
   Ne(0,15)+=reg20+reg32; Ne(1,16)+=reg20+reg32; Ne(2,17)+=reg20+reg32; Ne(0,6)+=reg39-reg36; Ne(1,7)+=reg39-reg36;
   Ne(2,8)+=reg39-reg36; Ne(0,0)+=reg46-reg24; Ne(1,1)+=reg46-reg24; Ne(2,2)+=reg46-reg24; Ne(0,12)+=reg38+reg20;
   Ne(1,13)+=reg38+reg20; Ne(2,14)+=reg38+reg20; Ne(0,9)+=reg45+reg15; Ne(1,10)+=reg45+reg15; Ne(2,11)+=reg45+reg15;
   Ne(0,3)+=reg21-reg24; Ne(1,4)+=reg21-reg24; Ne(2,5)+=reg21-reg24;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.12200846792814621817*e.pos(1)[0]; double reg1=0.36602540378443865451*e.pos(0)[0]; double reg2=0.36602540378443865451*e.pos(1)[0]; double reg3=0.12200846792814621817*e.pos(0)[0]; double reg4=0.12200846792814621817*e.pos(1)[1];
   double reg5=0.36602540378443865451*e.pos(0)[1]; double reg6=0.36602540378443865451*e.pos(1)[1]; double reg7=0.12200846792814621817*e.pos(0)[1]; double reg8=0.45534180126147951289*e.pos(0)[0]; double reg9=0.45534180126147951289*e.pos(0)[1];
   double reg10=0.12200846792814621817*e.pos(0)[2]; double reg11=reg0+reg1; double reg12=0.36602540378443865451*e.pos(1)[2]; double reg13=0.45534180126147951289*e.pos(2)[0]; double reg14=1.3660254037844385386*e.pos(1)[0];
   double reg15=1.3660254037844385386*e.pos(1)[1]; double reg16=reg4+reg5; double reg17=1.3660254037844385386*e.pos(2)[1]; double reg18=0.12200846792814621817*e.pos(1)[2]; double reg19=0.36602540378443865451*e.pos(0)[2];
   double reg20=0.45534180126147951289*e.pos(2)[1]; double reg21=reg6+reg7; double reg22=reg2+reg3; double reg23=1.3660254037844385386*e.pos(2)[0]; double reg24=reg13+reg22;
   double reg25=reg18+reg19; double reg26=1.3660254037844385386*e.pos(2)[2]; double reg27=0.12200846792814621817*e.pos(2)[0]; double reg28=1.3660254037844385386*e.pos(1)[2]; double reg29=reg21+reg20;
   double reg30=reg9+reg15; double reg31=0.45534180126147951289*e.pos(1)[0]; double reg32=1.3660254037844385386*e.pos(3)[1]; double reg33=reg8+reg14; double reg34=0.45534180126147951289*e.pos(1)[1];
   double reg35=1.3660254037844385386*e.pos(0)[1]; double reg36=reg12+reg10; double reg37=0.45534180126147951289*e.pos(2)[2]; double reg38=0.45534180126147951289*e.pos(0)[2]; double reg39=0.12200846792814621817*e.pos(2)[1];
   double reg40=1.3660254037844385386*e.pos(3)[0]; reg16=reg16+reg17; double reg41=0.45534180126147951289*e.pos(3)[0]; reg11=reg23+reg11; double reg42=1.3660254037844385386*e.pos(0)[0];
   double reg43=0.45534180126147951289*e.pos(3)[1]; double reg44=reg34+reg35; double reg45=reg38+reg28; double reg46=reg24+reg40; double reg47=0.48803387171258487271*e.pos(4)[0];
   double reg48=0.36602540378443865451*e.pos(3)[1]; reg30=reg39+reg30; double reg49=0.36602540378443865451*e.pos(2)[0]; double reg50=reg31+reg42; double reg51=reg29+reg32;
   reg33=reg27+reg33; double reg52=0.36602540378443865451*e.pos(3)[0]; double reg53=0.12200846792814621817*e.pos(2)[2]; double reg54=1.3660254037844385386*e.pos(3)[2]; double reg55=reg36+reg37;
   double reg56=0.48803387171258487271*e.pos(4)[1]; reg25=reg25+reg26; double reg57=0.45534180126147951289*e.pos(3)[2]; double reg58=0.45534180126147951289*e.pos(1)[2]; double reg59=1.3660254037844385386*e.pos(0)[2];
   reg11=reg11+reg41; reg16=reg16+reg43; double reg60=0.36602540378443865451*e.pos(2)[1]; reg11=reg11-reg47; double reg61=0.12200846792814621817*e.pos(3)[1];
   double reg62=reg49+reg50; double reg63=reg58+reg59; double reg64=0.36602540378443865451*e.pos(2)[2]; double reg65=0.66666666666666670528*e.pos(5)[0]; double reg66=reg55+reg54;
   double reg67=reg44+reg60; reg33=reg33+reg52; reg16=reg16-reg56; double reg68=0.66666666666666670528*e.pos(5)[1]; double reg69=reg47-reg46;
   reg25=reg25+reg57; double reg70=0.48803387171258487271*e.pos(4)[2]; double reg71=0.36602540378443865451*e.pos(3)[2]; reg45=reg53+reg45; double reg72=1.8213672050459180516*e.pos(4)[0];
   double reg73=reg56-reg51; double reg74=0.12200846792814621817*e.pos(3)[0]; reg30=reg30+reg48; double reg75=1.8213672050459180516*e.pos(4)[1]; reg33=reg33-reg72;
   reg16=reg16+reg68; double reg76=1.8213672050459180516*e.pos(6)[1]; double reg77=reg70-reg66; reg73=reg68+reg73; reg69=reg65+reg69;
   double reg78=1.8213672050459180516*e.pos(4)[2]; reg45=reg45+reg71; double reg79=reg63+reg64; double reg80=0.12200846792814621817*e.pos(3)[2]; double reg81=reg74+reg62;
   reg30=reg30-reg75; double reg82=1.8213672050459180516*e.pos(6)[0]; reg25=reg25-reg70; double reg83=0.66666666666666670528*e.pos(5)[2]; double reg84=reg61+reg67;
   reg11=reg65+reg11; reg69=reg82+reg69; reg77=reg83+reg77; reg73=reg76+reg73; reg11=reg11-reg82;
   double reg85=0.66666666666666670528*e.pos(7)[0]; reg16=reg16-reg76; double reg86=0.66666666666666670528*e.pos(7)[1]; reg45=reg45-reg78; reg25=reg25+reg83;
   double reg87=1.8213672050459180516*e.pos(6)[2]; double reg88=0.48803387171258487271*e.pos(6)[1]; reg30=reg68+reg30; double reg89=reg1+reg31; double reg90=0.48803387171258487271*e.pos(6)[0];
   reg33=reg65+reg33; double reg91=reg5+reg34; double reg92=reg72-reg81; double reg93=reg2+reg8; double reg94=reg75-reg84;
   double reg95=reg80+reg79; double reg96=reg6+reg9; reg15=reg7+reg15; reg14=reg3+reg14; reg27=reg27+reg93;
   reg77=reg87+reg77; reg3=0.48803387171258487271*e.pos(6)[2]; reg45=reg83+reg45; reg39=reg39+reg96; reg7=reg12+reg38;
   reg69=reg69-reg85; double reg97=reg78-reg95; reg30=reg30-reg88; reg94=reg68+reg94; reg92=reg65+reg92;
   reg33=reg33-reg90; reg73=reg73-reg86; reg65=reg19+reg58; reg11=reg11-reg85; reg17=reg17+reg91;
   reg16=reg16-reg86; reg68=0.66666666666666670528*e.pos(7)[2]; reg25=reg25-reg87; reg23=reg23+reg89; reg14=reg13+reg14;
   reg30=reg30-reg86; reg53=reg53+reg7; reg97=reg83+reg97; reg33=reg33-reg85; reg77=reg77-reg68;
   reg74=reg23+reg74; reg23=0.66666666666666670528*e.pos(4)[0]; reg27=reg40+reg27; reg61=reg17+reg61; reg17=0.66666666666666670528*e.pos(4)[1];
   reg15=reg20+reg15; reg92=reg90+reg92; reg26=reg26+reg65; reg45=reg45-reg3; reg42=reg0+reg42;
   reg0=pow(reg69,2); reg40=pow(reg11,2); reg83=pow(reg16,2); reg35=reg4+reg35; reg39=reg32+reg39;
   reg25=reg25-reg68; reg4=pow(reg73,2); reg94=reg88+reg94; reg28=reg10+reg28; reg10=1.8213672050459180516*e.pos(5)[1];
   reg32=reg52+reg14; reg39=reg39-reg17; reg80=reg26+reg80; reg85=reg92-reg85; reg26=0.48803387171258487271*e.pos(5)[0];
   reg27=reg27-reg23; reg86=reg94-reg86; reg92=0.48803387171258487271*e.pos(5)[1]; reg94=0.66666666666666670528*e.pos(4)[2]; double reg98=reg48+reg15;
   reg45=reg45-reg68; reg42=reg49+reg42; reg83=reg40+reg83; reg35=reg60+reg35; reg40=pow(reg25,2);
   reg53=reg54+reg53; reg4=reg0+reg4; reg28=reg37+reg28; reg0=pow(reg30,2); reg59=reg18+reg59;
   reg18=pow(reg33,2); reg54=pow(reg77,2); reg74=reg74-reg23; double reg99=1.8213672050459180516*e.pos(5)[0]; reg97=reg97+reg3;
   reg61=reg61-reg17; double reg100=reg23+reg32; reg80=reg80-reg94; double reg101=0.48803387171258487271*e.pos(5)[2]; double reg102=pow(reg86,2);
   double reg103=0.66666666666666670528*e.pos(6)[1]; reg61=reg61-reg10; double reg104=pow(reg85,2); double reg105=0.66666666666666670528*e.pos(6)[0]; reg74=reg74-reg99;
   reg40=reg83+reg40; reg0=reg18+reg0; reg18=pow(reg45,2); reg68=reg97-reg68; reg83=reg17+reg98;
   reg27=reg27-reg26; reg54=reg4+reg54; reg59=reg64+reg59; reg39=reg39-reg92; reg4=reg43+reg35;
   reg97=reg71+reg28; reg53=reg53-reg94; double reg106=reg41+reg42; double reg107=1.8213672050459180516*e.pos(5)[2]; double reg108=pow(reg68,2);
   reg53=reg53-reg101; reg27=reg105+reg27; double reg109=1.8213672050459180516*e.pos(7)[0]; double reg110=reg57+reg59; double reg111=reg10-reg83;
   reg61=reg61+reg103; double reg112=0.48803387171258487271*e.pos(7)[1]; reg80=reg80-reg107; reg18=reg0+reg18; reg0=0.48803387171258487271*e.pos(7)[0];
   reg40=pow(reg40,0.5); reg23=reg23+reg106; reg102=reg104+reg102; reg74=reg74+reg105; reg54=pow(reg54,0.5);
   reg39=reg103+reg39; reg104=1.8213672050459180516*e.pos(7)[1]; reg17=reg17+reg4; double reg113=reg94+reg97; double reg114=reg99-reg100;
   double reg115=0.66666666666666670528*e.pos(6)[2]; double reg116=reg69/reg54; reg111=reg103+reg111; reg94=reg94+reg110; reg18=pow(reg18,0.5);
   double reg117=reg16/reg40; double reg118=reg11/reg40; double reg119=reg92-reg17; double reg120=reg26-reg23; double reg121=reg107-reg113;
   reg80=reg80+reg115; double reg122=0.48803387171258487271*e.pos(7)[2]; reg39=reg39-reg104; reg61=reg61-reg112; reg114=reg105+reg114;
   reg27=reg27-reg109; double reg123=1.8213672050459180516*e.pos(7)[2]; reg53=reg115+reg53; reg74=reg74-reg0; reg108=reg102+reg108;
   reg102=reg73/reg54; reg53=reg53-reg123; reg114=reg0+reg114; reg120=reg105+reg120; reg105=reg102*reg39;
   reg119=reg103+reg119; reg121=reg115+reg121; reg111=reg112+reg111; reg80=reg80-reg122; reg103=reg117*reg61;
   double reg124=reg116*reg27; double reg125=reg33/reg18; reg40=reg25/reg40; double reg126=reg30/reg18; double reg127=reg101-reg94;
   reg54=reg77/reg54; double reg128=reg118*reg74; reg108=pow(reg108,0.5); double reg129=reg126*reg111; reg18=reg45/reg18;
   double reg130=reg125*reg114; reg121=reg122+reg121; double reg131=reg40*reg80; reg105=reg124+reg105; reg124=reg54*reg53;
   reg127=reg115+reg127; reg115=reg85/reg108; reg103=reg128+reg103; reg128=reg86/reg108; reg119=reg104+reg119;
   reg120=reg109+reg120; reg131=reg103+reg131; reg124=reg105+reg124; reg103=reg128*reg119; reg105=reg115*reg120;
   reg108=reg68/reg108; double reg132=reg18*reg121; reg127=reg123+reg127; reg129=reg130+reg129; reg130=reg118*reg131;
   double reg133=reg117*reg131; reg103=reg105+reg103; reg105=reg116*reg124; double reg134=reg102*reg124; double reg135=reg108*reg127;
   reg132=reg129+reg132; reg129=reg125*reg132; double reg136=reg40*reg131; reg133=reg61-reg133; reg130=reg74-reg130;
   reg105=reg27-reg105; reg134=reg39-reg134; double reg137=reg54*reg124; double reg138=reg126*reg132; reg135=reg103+reg135;
   reg138=reg111-reg138; reg103=reg18*reg132; double reg139=reg128*reg135; double reg140=reg115*reg135; reg137=reg53-reg137;
   double reg141=pow(reg134,2); double reg142=pow(reg105,2); reg136=reg80-reg136; double reg143=pow(reg133,2); double reg144=pow(reg130,2);
   reg129=reg114-reg129; reg103=reg121-reg103; double reg145=pow(reg136,2); reg143=reg144+reg143; reg144=pow(reg137,2);
   double reg146=pow(reg138,2); double reg147=reg108*reg135; reg139=reg119-reg139; double reg148=pow(reg129,2); reg140=reg120-reg140;
   reg141=reg142+reg141; reg145=reg143+reg145; reg146=reg148+reg146; reg142=pow(reg140,2); reg143=pow(reg139,2);
   reg147=reg127-reg147; reg144=reg141+reg144; reg141=pow(reg103,2); reg145=pow(reg145,0.5); reg144=pow(reg144,0.5);
   reg148=pow(reg147,2); reg141=reg146+reg141; reg143=reg142+reg143; reg134=reg134/reg144; reg105=reg105/reg144;
   reg141=pow(reg141,0.5); reg133=reg133/reg145; reg130=reg130/reg145; reg148=reg143+reg148; reg142=reg11*reg130;
   reg143=reg16*reg133; reg116=reg69*reg116; reg102=reg73*reg102; reg148=pow(reg148,0.5); reg117=reg16*reg117;
   reg118=reg11*reg118; reg145=reg136/reg145; reg133=reg61*reg133; reg130=reg74*reg130; reg27=reg27*reg105;
   reg138=reg138/reg141; reg129=reg129/reg141; reg73=reg73*reg134; reg105=reg69*reg105; reg144=reg137/reg144;
   reg134=reg39*reg134; reg40=reg25*reg40; reg117=reg118+reg117; reg143=reg142+reg143; reg25=reg25*reg145;
   reg141=reg103/reg141; reg145=reg80*reg145; reg111=reg111*reg138; reg114=reg114*reg129; reg134=reg27+reg134;
   reg129=reg33*reg129; reg133=reg130+reg133; reg140=reg140/reg148; reg125=reg33*reg125; reg139=reg139/reg148;
   reg138=reg30*reg138; reg102=reg116+reg102; reg54=reg77*reg54; reg126=reg30*reg126; reg77=reg77*reg144;
   reg73=reg105+reg73; reg144=reg53*reg144; reg40=reg117+reg40; reg138=reg129+reg138; reg11=reg45*reg141;
   reg141=reg121*reg141; reg111=reg114+reg111; reg126=reg125+reg126; reg18=reg45*reg18; reg120=reg120*reg140;
   reg119=reg119*reg139; reg148=reg147/reg148; reg140=reg85*reg140; reg139=reg86*reg139; reg77=reg73+reg77;
   reg115=reg85*reg115; reg128=reg86*reg128; reg145=reg133+reg145; reg54=reg102+reg54; reg144=reg134+reg144;
   reg25=reg143+reg25; reg77=reg124*reg77; reg18=reg126+reg18; reg119=reg120+reg119; reg128=reg115+reg128;
   reg108=reg68*reg108; reg127=reg127*reg148; reg25=reg131*reg25; reg144=reg54*reg144; reg141=reg111+reg141;
   reg139=reg140+reg139; reg145=reg40*reg145; reg148=reg68*reg148; reg11=reg138+reg11; reg127=reg119+reg127;
   reg148=reg139+reg148; reg25=reg145-reg25; reg108=reg128+reg108; reg77=reg144-reg77; reg141=reg18*reg141;
   reg11=reg132*reg11; reg16=0.024056261216234395431*reg77; reg18=0.024056261216234395431*reg25; reg27=0.035220810900864524453*reg25; reg30=0.035220810900864524453*reg77;
   reg33=0.024056261216234409915*reg77; reg11=reg141-reg11; reg39=0.04166666666666666908*reg25; reg40=0.13144585576580215187*reg25; reg148=reg135*reg148;
   reg45=0.13144585576580215187*reg77; reg53=0.04166666666666666908*reg77; reg127=reg108*reg127; reg54=0.024056261216234409915*reg25; reg61=0.13144585576580215187*reg11;
   reg68=reg40+reg45; reg40=reg30+reg40; reg30=reg27+reg30; reg69=0.035220810900864524453*reg11; reg73=0.024056261216234409915*reg11;
   reg16=reg16-reg39; reg74=0.04166666666666666908*reg11; reg54=reg54+reg53; reg148=reg127-reg148; reg80=0.024056261216234395431*reg11;
   reg33=reg39+reg33; reg45=reg27+reg45; reg53=reg18-reg53; reg54=reg54+reg74; reg18=0.024056261216234409915*reg148;
   reg27=0.035220810900864524453*reg148; reg40=reg61+reg40; reg74=reg53-reg74; reg39=0.024056261216234395431*reg148; reg68=reg68+reg69;
   reg53=0.13144585576580215187*reg148; reg61=reg30+reg61; reg73=reg16-reg73; reg33=reg80-reg33; reg16=0.04166666666666666908*reg148;
   reg45=reg69+reg45; Ne(0,21)+=reg53+reg45; Ne(1,22)+=reg53+reg45; Ne(2,23)+=reg53+reg45; Ne(0,18)+=reg27+reg68;
   Ne(1,19)+=reg27+reg68; Ne(2,20)+=reg27+reg68; Ne(0,15)+=reg40+reg27; Ne(1,16)+=reg40+reg27; Ne(2,17)+=reg40+reg27;
   Ne(0,12)+=reg61+reg53; Ne(1,13)+=reg61+reg53; Ne(2,14)+=reg61+reg53; Ne(0,9)+=reg73-reg16; Ne(1,10)+=reg73-reg16;
   Ne(2,11)+=reg73-reg16; Ne(0,6)+=reg74-reg18; Ne(1,7)+=reg74-reg18; Ne(2,8)+=reg74-reg18; Ne(0,3)+=reg33-reg16;
   Ne(1,4)+=reg33-reg16; Ne(2,5)+=reg33-reg16; Ne(0,0)+=reg39-reg54; Ne(1,1)+=reg39-reg54; Ne(2,2)+=reg39-reg54;

}

};
} // namespace LMT

