
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
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.56758792732771999991*e.pos(1)[0]; double reg3=0.78379396366385999995*e.pos(0)[0]; double reg4=0.78379396366385999995*e.pos(1)[1];
   double reg5=0.56758792732771999991*e.pos(0)[1]; double reg6=0.78379396366385999995*e.pos(0)[1]; double reg7=0.56758792732771999991*e.pos(1)[1]; double reg8=reg4+reg5; double reg9=1.3513818909915799999*e.pos(3)[1];
   double reg10=reg2+reg3; double reg11=0.78379396366385999995*e.pos(1)[2]; double reg12=1.3513818909915799999*e.pos(3)[0]; double reg13=reg0+reg1; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg7+reg6; double reg16=0.56758792732771999991*e.pos(1)[2]; double reg17=0.78379396366385999995*e.pos(0)[2]; double reg18=reg11+reg14; double reg19=0.63369514596091600003*e.pos(1)[0];
   double reg20=1.3513818909915799999*e.pos(3)[2]; double reg21=reg12-reg10; double reg22=2.2673902919218320001*e.pos(0)[1]; double reg23=0.63369514596091600003*e.pos(1)[1]; double reg24=1.78379396366386*e.pos(4)[1];
   reg8=reg8-reg9; double reg25=2.2673902919218320001*e.pos(0)[0]; double reg26=reg9-reg15; double reg27=1.78379396366386*e.pos(4)[0]; reg13=reg13-reg12;
   double reg28=reg16+reg17; double reg29=reg23+reg22; reg21=reg27+reg21; double reg30=0.63369514596091600003*e.pos(1)[2]; double reg31=2.2673902919218320001*e.pos(0)[2];
   double reg32=2.9010854378827480001*e.pos(3)[1]; double reg33=0.63369514596091600003*e.pos(0)[1]; double reg34=2.2673902919218320001*e.pos(1)[1]; double reg35=0.63369514596091600003*e.pos(0)[0]; reg4=reg4-reg6;
   reg0=reg0-reg3; double reg36=0.43241207267228000009*e.pos(4)[1]; double reg37=2.2673902919218320001*e.pos(1)[0]; double reg38=reg20-reg28; double reg39=0.43241207267228000009*e.pos(4)[0];
   double reg40=reg19+reg25; reg13=reg13+reg27; double reg41=1.78379396366386*e.pos(5)[0]; reg8=reg8+reg24; double reg42=1.78379396366386*e.pos(5)[1];
   double reg43=2.9010854378827480001*e.pos(3)[0]; reg26=reg24+reg26; reg18=reg18-reg20; double reg44=1.78379396366386*e.pos(4)[2]; double reg45=reg32-reg29;
   double reg46=reg43-reg40; double reg47=0.366304854039084*e.pos(4)[1]; reg4=reg36+reg4; double reg48=0.366304854039084*e.pos(4)[0]; double reg49=2.9010854378827480001*e.pos(3)[2];
   reg11=reg11-reg17; double reg50=0.43241207267228000009*e.pos(5)[1]; reg38=reg44+reg38; reg13=reg13-reg41; reg26=reg26-reg42;
   reg42=reg8-reg42; reg18=reg18+reg44; reg8=1.78379396366386*e.pos(5)[2]; reg41=reg21-reg41; reg21=0.63369514596091600003*e.pos(0)[2];
   double reg51=2.2673902919218320001*e.pos(1)[2]; reg34=reg34+reg33; reg37=reg37+reg35; reg0=reg39+reg0; double reg52=0.43241207267228000009*e.pos(4)[2];
   double reg53=0.43241207267228000009*e.pos(5)[0]; double reg54=reg30+reg31; double reg55=pow(reg41,2); reg50=reg4-reg50; reg11=reg52+reg11;
   reg4=pow(reg26,2); reg53=reg0-reg53; reg38=reg38-reg8; reg0=0.43241207267228000009*e.pos(5)[2]; reg34=reg34-reg32;
   double reg56=0.78379396366385999995*e.pos(2)[0]; double reg57=pow(reg42,2); reg37=reg37-reg43; double reg58=pow(reg13,2); double reg59=2.710505431213761085e-20*e.pos(3)[0];
   double reg60=0.78379396366385999995*e.pos(2)[1]; double reg61=0.366304854039084*e.pos(4)[2]; double reg62=reg49-reg54; reg8=reg18-reg8; reg18=reg33-reg23;
   double reg63=0.366304854039084*e.pos(5)[1]; reg45=reg45+reg47; double reg64=2.710505431213761085e-20*e.pos(3)[1]; double reg65=0.366304854039084*e.pos(5)[0]; reg51=reg51+reg21;
   reg46=reg46+reg48; double reg66=reg35-reg19; double reg67=3.2673902919218320001*e.pos(4)[1]; reg66=reg66-reg59; reg18=reg18-reg64;
   double reg68=2.710505431213761085e-20*e.pos(3)[2]; double reg69=reg21-reg30; double reg70=3.2673902919218320001*e.pos(4)[0]; reg51=reg51-reg49; reg34=reg47+reg34;
   reg37=reg48+reg37; double reg71=0.366304854039084*e.pos(5)[2]; reg62=reg62+reg61; reg45=reg45-reg63; double reg72=pow(reg53,2);
   reg46=reg46-reg65; double reg73=pow(reg50,2); double reg74=0.56758792732771999991*e.pos(2)[1]; reg0=reg11-reg0; reg11=0.56758792732771999991*e.pos(2)[0];
   double reg75=reg1+reg56; double reg76=1.78379396366386*e.pos(3)[0]; reg4=reg55+reg4; reg55=pow(reg8,2); reg57=reg58+reg57;
   reg58=reg5+reg60; double reg77=1.78379396366386*e.pos(3)[1]; double reg78=0.78379396366385999995*e.pos(2)[2]; double reg79=pow(reg38,2); double reg80=3.2673902919218320001*e.pos(5)[0];
   reg18=reg67+reg18; reg66=reg70+reg66; double reg81=3.2673902919218320001*e.pos(5)[1]; double reg82=3.2673902919218320001*e.pos(4)[2]; reg55=reg57+reg55;
   reg69=reg69-reg68; reg73=reg72+reg73; reg51=reg61+reg51; reg57=pow(reg0,2); reg63=reg34-reg63;
   reg34=reg3+reg11; reg65=reg37-reg65; reg37=reg6+reg74; reg75=reg75-reg76; reg72=0.56758792732771999991*e.pos(2)[2];
   double reg83=0.63369514596091600003*e.pos(2)[1]; double reg84=0.63369514596091600003*e.pos(2)[0]; reg62=reg62-reg71; reg58=reg58-reg77; double reg85=pow(reg45,2);
   double reg86=reg14+reg78; double reg87=1.78379396366386*e.pos(3)[2]; double reg88=pow(reg46,2); reg79=reg4+reg79; reg4=0.43241207267228000009*e.pos(3)[1];
   reg6=reg60-reg6; reg3=reg56-reg3; reg56=0.43241207267228000009*e.pos(3)[0]; reg71=reg51-reg71; reg57=reg73+reg57;
   reg51=pow(reg63,2); reg4=reg6-reg4; reg76=reg76+reg34; reg6=pow(reg65,2); reg77=reg77+reg37;
   reg56=reg3-reg56; reg3=reg17+reg72; reg75=reg27+reg75; reg60=1.3513818909915799999*e.pos(5)[0]; reg73=0.63369514596091600003*e.pos(2)[2];
   double reg89=0.366304854039084*e.pos(3)[1]; reg22=reg22+reg83; double reg90=0.366304854039084*e.pos(3)[0]; reg25=reg25+reg84; double reg91=pow(reg62,2);
   reg58=reg24+reg58; double reg92=1.3513818909915799999*e.pos(5)[1]; reg85=reg88+reg85; reg79=pow(reg79,0.5); reg86=reg86-reg87;
   reg17=reg78-reg17; reg78=0.43241207267228000009*e.pos(3)[2]; reg80=reg66-reg80; reg81=reg18-reg81; reg55=pow(reg55,0.5);
   reg69=reg82+reg69; reg18=3.2673902919218320001*e.pos(5)[2]; reg24=reg24-reg77; reg39=reg56+reg39; reg87=reg87+reg3;
   reg56=1.3513818909915799999*e.pos(5)[2]; reg66=0.366304854039084*e.pos(3)[2]; reg31=reg31+reg73; reg78=reg17-reg78; reg17=reg22+reg89;
   reg88=2.2673902919218320001*e.pos(2)[0]; double reg93=reg25+reg90; reg86=reg44+reg86; double reg94=pow(reg80,2); double reg95=2.2673902919218320001*e.pos(2)[1];
   reg91=reg85+reg91; reg75=reg75-reg60; reg18=reg69-reg18; reg58=reg58-reg92; reg69=reg41/reg79;
   reg85=pow(reg81,2); double reg96=reg26/reg79; reg57=pow(reg57,0.5); double reg97=reg33-reg83; reg36=reg4+reg36;
   reg4=3.2673902919218320001*e.pos(3)[1]; double reg98=pow(reg71,2); reg27=reg27-reg76; double reg99=reg13/reg55; double reg100=reg42/reg55;
   double reg101=reg35-reg84; double reg102=3.2673902919218320001*e.pos(3)[0]; reg51=reg6+reg51; reg91=pow(reg91,0.5); reg6=reg50/reg57;
   double reg103=reg53/reg57; double reg104=pow(reg18,2); reg95=reg33+reg95; double reg105=reg48-reg93; double reg106=reg99*reg75;
   reg79=reg38/reg79; double reg107=2.2673902919218320001*e.pos(2)[2]; double reg108=reg96*reg36; reg98=reg51+reg98; reg51=reg100*reg58;
   reg88=reg35+reg88; reg85=reg94+reg85; reg86=reg86-reg56; reg102=reg101-reg102; reg94=reg69*reg39;
   reg4=reg97-reg4; reg24=reg92+reg24; reg27=reg60+reg27; reg97=reg21-reg73; reg101=3.2673902919218320001*e.pos(3)[2];
   reg44=reg44-reg87; double reg109=reg31+reg66; reg55=reg8/reg55; double reg110=2.9010854378827480001*e.pos(5)[0]; reg52=reg78+reg52;
   reg78=2.9010854378827480001*e.pos(5)[1]; double reg111=reg47-reg17; double reg112=5.42101086242752217e-20*e.pos(5)[0]; reg67=reg4+reg67; reg105=reg105+reg110;
   reg104=reg85+reg104; reg98=pow(reg98,0.5); reg111=reg111+reg78; reg4=5.42101086242752217e-20*e.pos(5)[1]; reg107=reg21+reg107;
   reg51=reg106+reg51; reg90=reg88-reg90; reg85=reg45/reg91; reg88=reg46/reg91; reg101=reg97-reg101;
   reg97=2.9010854378827480001*e.pos(5)[2]; reg70=reg102+reg70; reg102=reg61-reg109; reg89=reg95-reg89; reg95=reg55*reg86;
   reg108=reg94+reg108; reg57=reg0/reg57; reg94=reg103*reg27; reg106=reg6*reg24; reg44=reg56+reg44;
   double reg113=reg79*reg52; reg66=reg107-reg66; reg91=reg62/reg91; reg89=reg47+reg89; reg90=reg48+reg90;
   reg47=reg88*reg105; reg48=reg85*reg111; reg102=reg102+reg97; reg104=pow(reg104,0.5); reg113=reg108+reg113;
   reg107=reg57*reg44; reg106=reg94+reg106; reg94=reg65/reg98; reg108=reg63/reg98; double reg114=5.42101086242752217e-20*e.pos(5)[2];
   reg82=reg101+reg82; reg67=reg67-reg4; reg70=reg70-reg112; reg95=reg51+reg95; reg51=reg69*reg113;
   reg101=reg80/reg104; double reg115=reg96*reg113; double reg116=reg81/reg104; reg48=reg47+reg48; reg47=reg100*reg95;
   reg90=reg90-reg110; double reg117=reg91*reg102; reg98=reg71/reg98; reg89=reg89-reg78; reg82=reg82-reg114;
   reg107=reg106+reg107; reg106=reg108*reg67; double reg118=reg99*reg95; double reg119=reg94*reg70; reg66=reg61+reg66;
   reg61=reg103*reg107; reg117=reg48+reg117; reg104=reg18/reg104; reg51=reg39-reg51; reg115=reg36-reg115;
   reg48=reg55*reg95; double reg120=reg79*reg113; reg47=reg58-reg47; double reg121=reg98*reg82; reg118=reg75-reg118;
   reg106=reg119+reg106; reg66=reg66-reg97; reg119=reg6*reg107; double reg122=reg116*reg89; double reg123=reg101*reg90;
   double reg124=reg104*reg66; reg122=reg123+reg122; reg123=pow(reg118,2); reg121=reg106+reg121; reg106=pow(reg47,2);
   reg120=reg52-reg120; double reg125=pow(reg115,2); reg48=reg86-reg48; double reg126=pow(reg51,2); double reg127=reg57*reg107;
   reg61=reg27-reg61; reg119=reg24-reg119; double reg128=reg88*reg117; double reg129=reg85*reg117; double reg130=pow(reg61,2);
   reg124=reg122+reg124; reg128=reg105-reg128; reg122=reg94*reg121; double reg131=reg108*reg121; reg126=reg125+reg126;
   reg106=reg123+reg106; reg129=reg111-reg129; reg123=pow(reg48,2); reg125=pow(reg120,2); double reg132=reg91*reg117;
   reg127=reg44-reg127; double reg133=pow(reg119,2); double reg134=pow(reg128,2); double reg135=reg98*reg121; reg131=reg67-reg131;
   double reg136=reg116*reg124; double reg137=pow(reg129,2); double reg138=pow(reg127,2); reg123=reg106+reg123; reg122=reg70-reg122;
   reg133=reg130+reg133; reg106=reg101*reg124; reg125=reg126+reg125; reg132=reg102-reg132; reg106=reg90-reg106;
   reg126=pow(reg122,2); reg136=reg89-reg136; reg130=pow(reg132,2); double reg139=pow(reg131,2); reg135=reg82-reg135;
   reg137=reg134+reg137; reg138=reg133+reg138; reg133=reg104*reg124; reg125=pow(reg125,0.5); reg123=pow(reg123,0.5);
   reg130=reg137+reg130; reg134=pow(reg136,2); reg133=reg66-reg133; reg51=reg51/reg125; reg115=reg115/reg125;
   reg137=pow(reg106,2); reg138=pow(reg138,0.5); reg47=reg47/reg123; reg118=reg118/reg123; reg139=reg126+reg139;
   reg126=pow(reg135,2); reg36=reg36*reg115; double reg140=pow(reg133,2); reg39=reg39*reg51; reg125=reg120/reg125;
   reg134=reg137+reg134; reg115=reg26*reg115; reg96=reg26*reg96; reg69=reg41*reg69; reg126=reg139+reg126;
   reg51=reg41*reg51; reg123=reg48/reg123; reg58=reg58*reg47; reg99=reg13*reg99; reg100=reg42*reg100;
   reg130=pow(reg130,0.5); reg61=reg61/reg138; reg75=reg75*reg118; reg47=reg42*reg47; reg119=reg119/reg138;
   reg118=reg13*reg118; reg58=reg75+reg58; reg103=reg53*reg103; reg6=reg50*reg6; reg126=pow(reg126,0.5);
   reg79=reg38*reg79; reg24=reg24*reg119; reg27=reg27*reg61; reg128=reg128/reg130; reg129=reg129/reg130;
   reg36=reg39+reg36; reg55=reg8*reg55; reg140=reg134+reg140; reg100=reg99+reg100; reg52=reg52*reg125;
   reg119=reg50*reg119; reg61=reg53*reg61; reg8=reg8*reg123; reg47=reg118+reg47; reg125=reg38*reg125;
   reg115=reg51+reg115; reg138=reg127/reg138; reg96=reg69+reg96; reg123=reg86*reg123; reg52=reg36+reg52;
   reg13=reg45*reg129; reg8=reg47+reg8; reg26=reg46*reg128; reg36=reg0*reg138; reg57=reg0*reg57;
   reg140=pow(reg140,0.5); reg119=reg61+reg119; reg130=reg132/reg130; reg129=reg111*reg129; reg123=reg58+reg123;
   reg128=reg105*reg128; reg88=reg46*reg88; reg125=reg115+reg125; reg138=reg44*reg138; reg24=reg27+reg24;
   reg79=reg96+reg79; reg131=reg131/reg126; reg55=reg100+reg55; reg122=reg122/reg126; reg85=reg45*reg85;
   reg6=reg103+reg6; reg13=reg26+reg13; reg0=reg62*reg130; reg130=reg102*reg130; reg129=reg128+reg129;
   reg91=reg62*reg91; reg85=reg88+reg85; reg138=reg24+reg138; reg36=reg119+reg36; reg126=reg135/reg126;
   reg67=reg67*reg131; reg70=reg70*reg122; reg94=reg65*reg94; reg122=reg65*reg122; reg131=reg63*reg131;
   reg108=reg63*reg108; reg123=reg55*reg123; reg136=reg136/reg140; reg125=reg113*reg125; reg106=reg106/reg140;
   reg52=reg79*reg52; reg57=reg6+reg57; reg8=reg95*reg8; reg125=reg52-reg125; reg138=reg57*reg138;
   reg98=reg71*reg98; reg91=reg85+reg91; reg108=reg94+reg108; reg8=reg123-reg8; reg36=reg107*reg36;
   reg67=reg70+reg67; reg0=reg13+reg0; reg6=reg81*reg136; reg13=reg80*reg106; reg130=reg129+reg130;
   reg106=reg90*reg106; reg136=reg89*reg136; reg101=reg80*reg101; reg71=reg71*reg126; reg116=reg81*reg116;
   reg131=reg122+reg131; reg140=reg133/reg140; reg126=reg82*reg126; reg66=reg66*reg140; reg136=reg106+reg136;
   reg24=0.009463616120767210603*reg125; reg26=0.088847818743090689935*reg8; reg27=0.088847818743090689935*reg125; reg38=0.005384432036113586778*reg8; reg140=reg18*reg140;
   reg39=0.021537728144454347112*reg125; reg41=0.021537728144454347112*reg8; reg6=reg13+reg6; reg104=reg18*reg104; reg116=reg101+reg116;
   reg13=0.009463616120767210603*reg8; reg0=reg117*reg0; reg130=reg91*reg130; reg36=reg138-reg36; reg71=reg131+reg71;
   reg18=0.005384432036113586778*reg125; reg126=reg67+reg126; reg98=reg108+reg98; reg27=reg41+reg27; reg24=reg38+reg24;
   reg38=reg18+reg38; reg42=0.009463616120767210603*reg36; reg44=0.021537728144454347112*reg36; reg26=reg39+reg26; reg45=0.005384432036113586778*reg36;
   reg39=reg41+reg39; reg18=reg13+reg18; reg13=0.088847818743090689935*reg36; reg104=reg116+reg104; reg140=reg6+reg140;
   reg0=reg130-reg0; reg71=reg121*reg71; reg126=reg98*reg126; reg66=reg136+reg66; reg6=0.016449618187943419918*reg0;
   reg13=reg39+reg13; reg71=reg126-reg71; reg18=reg18+reg45; reg66=reg104*reg66; reg42=reg38+reg42;
   reg26=reg26+reg44; reg38=0.028457289286966203713*reg0; reg39=0.0018441552587796664112*reg0; reg41=0.0041124045469858549794*reg0; reg24=reg45+reg24;
   reg140=reg124*reg140; reg27=reg44+reg27; reg27=reg6+reg27; reg39=reg26+reg39; reg26=0.0018441552587796664109*reg71;
   reg44=0.016449618187943419918*reg71; reg45=0.016449618187943419916*reg71; reg46=0.004112404546985854979*reg71; reg18=reg38-reg18; reg140=reg66-reg140;
   reg38=0.028457289286966203713*reg71; reg24=reg24+reg41; reg42=reg41+reg42; reg41=0.0041124045469858549794*reg71; reg6=reg13+reg6;
   reg13=0.028457289286966203713*reg140; reg41=reg42+reg41; reg24=reg38-reg24; reg26=reg27+reg26; reg27=0.0018441552587796664111*reg140;
   reg44=reg39+reg44; reg38=0.016449618187943419918*reg140; reg46=reg18-reg46; reg18=0.0041124045469858549794*reg140; reg45=reg6+reg45;
   Ne(0,15)+=reg38+reg26; Ne(1,16)+=reg38+reg26; Ne(2,17)+=reg38+reg26; Ne(0,6)+=reg13-reg41; Ne(1,7)+=reg13-reg41;
   Ne(2,8)+=reg13-reg41; Ne(0,0)+=reg46-reg18; Ne(1,1)+=reg46-reg18; Ne(2,2)+=reg46-reg18; Ne(0,12)+=reg44+reg38;
   Ne(1,13)+=reg44+reg38; Ne(2,14)+=reg44+reg38; Ne(0,9)+=reg45+reg27; Ne(1,10)+=reg45+reg27; Ne(2,11)+=reg45+reg27;
   Ne(0,3)+=reg24-reg18; Ne(1,4)+=reg24-reg18; Ne(2,5)+=reg24-reg18;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(0)[0]; double reg1=0.36602540378443865451*e.pos(1)[0]; double reg2=0.12200846792814621817*e.pos(0)[0]; double reg3=0.12200846792814621817*e.pos(1)[1]; double reg4=0.36602540378443865451*e.pos(0)[1];
   double reg5=0.12200846792814621817*e.pos(1)[0]; double reg6=0.36602540378443865451*e.pos(1)[1]; double reg7=0.12200846792814621817*e.pos(0)[1]; double reg8=0.45534180126147951289*e.pos(0)[0]; double reg9=0.45534180126147951289*e.pos(0)[1];
   double reg10=0.12200846792814621817*e.pos(0)[2]; double reg11=0.36602540378443865451*e.pos(1)[2]; double reg12=1.3660254037844385386*e.pos(1)[0]; double reg13=reg5+reg0; double reg14=1.3660254037844385386*e.pos(1)[1];
   double reg15=0.45534180126147951289*e.pos(2)[1]; double reg16=reg6+reg7; double reg17=0.45534180126147951289*e.pos(2)[0]; double reg18=0.36602540378443865451*e.pos(0)[2]; double reg19=0.12200846792814621817*e.pos(1)[2];
   double reg20=1.3660254037844385386*e.pos(2)[1]; double reg21=reg3+reg4; double reg22=1.3660254037844385386*e.pos(2)[0]; double reg23=reg1+reg2; double reg24=reg19+reg18;
   double reg25=1.3660254037844385386*e.pos(2)[2]; double reg26=1.3660254037844385386*e.pos(1)[2]; double reg27=reg16+reg15; double reg28=reg9+reg14; double reg29=0.45534180126147951289*e.pos(1)[0];
   double reg30=1.3660254037844385386*e.pos(3)[1]; double reg31=reg8+reg12; double reg32=0.45534180126147951289*e.pos(1)[1]; double reg33=1.3660254037844385386*e.pos(0)[1]; double reg34=reg11+reg10;
   double reg35=0.45534180126147951289*e.pos(2)[2]; double reg36=0.45534180126147951289*e.pos(0)[2]; double reg37=0.12200846792814621817*e.pos(2)[1]; double reg38=0.12200846792814621817*e.pos(2)[0]; double reg39=reg17+reg23;
   double reg40=1.3660254037844385386*e.pos(0)[0]; reg13=reg22+reg13; double reg41=0.45534180126147951289*e.pos(3)[0]; double reg42=1.3660254037844385386*e.pos(3)[0]; reg21=reg21+reg20;
   double reg43=0.45534180126147951289*e.pos(3)[1]; reg13=reg13+reg41; double reg44=0.36602540378443865451*e.pos(2)[1]; double reg45=reg36+reg26; double reg46=reg39+reg42;
   double reg47=reg32+reg33; double reg48=0.48803387171258487271*e.pos(4)[0]; double reg49=reg27+reg30; double reg50=0.36602540378443865451*e.pos(3)[1]; reg28=reg37+reg28;
   double reg51=0.36602540378443865451*e.pos(2)[0]; double reg52=1.3660254037844385386*e.pos(0)[2]; double reg53=0.45534180126147951289*e.pos(1)[2]; double reg54=reg29+reg40; reg31=reg38+reg31;
   double reg55=0.45534180126147951289*e.pos(3)[2]; reg24=reg24+reg25; reg21=reg21+reg43; double reg56=reg34+reg35; double reg57=1.3660254037844385386*e.pos(3)[2];
   double reg58=0.36602540378443865451*e.pos(3)[0]; double reg59=0.12200846792814621817*e.pos(2)[2]; double reg60=0.48803387171258487271*e.pos(4)[1]; reg31=reg31+reg58; reg13=reg13-reg48;
   double reg61=reg51+reg54; double reg62=0.12200846792814621817*e.pos(3)[1]; double reg63=reg47+reg44; double reg64=reg53+reg52; double reg65=0.36602540378443865451*e.pos(2)[2];
   double reg66=reg56+reg57; double reg67=0.66666666666666670528*e.pos(5)[0]; reg21=reg21-reg60; double reg68=0.66666666666666670528*e.pos(5)[1]; double reg69=reg48-reg46;
   reg24=reg24+reg55; double reg70=0.48803387171258487271*e.pos(4)[2]; double reg71=0.36602540378443865451*e.pos(3)[2]; reg45=reg59+reg45; double reg72=1.8213672050459180516*e.pos(4)[1];
   double reg73=1.8213672050459180516*e.pos(4)[0]; double reg74=reg60-reg49; double reg75=0.12200846792814621817*e.pos(3)[0]; reg28=reg28+reg50; reg31=reg31-reg73;
   reg21=reg21+reg68; double reg76=1.8213672050459180516*e.pos(6)[1]; double reg77=reg70-reg66; reg74=reg68+reg74; reg69=reg67+reg69;
   double reg78=1.8213672050459180516*e.pos(4)[2]; reg45=reg45+reg71; double reg79=reg64+reg65; double reg80=0.12200846792814621817*e.pos(3)[2]; reg13=reg67+reg13;
   double reg81=reg75+reg61; reg24=reg24-reg70; double reg82=0.66666666666666670528*e.pos(5)[2]; double reg83=reg62+reg63; double reg84=1.8213672050459180516*e.pos(6)[0];
   reg28=reg28-reg72; reg69=reg84+reg69; reg77=reg82+reg77; reg74=reg76+reg74; reg13=reg13-reg84;
   double reg85=0.66666666666666670528*e.pos(7)[0]; reg21=reg21-reg76; double reg86=0.66666666666666670528*e.pos(7)[1]; reg45=reg45-reg78; reg24=reg24+reg82;
   double reg87=1.8213672050459180516*e.pos(6)[2]; double reg88=0.48803387171258487271*e.pos(6)[1]; reg28=reg68+reg28; double reg89=reg0+reg29; double reg90=0.48803387171258487271*e.pos(6)[0];
   reg31=reg67+reg31; double reg91=reg4+reg32; double reg92=reg73-reg81; double reg93=reg1+reg8; double reg94=reg72-reg83;
   double reg95=reg80+reg79; double reg96=reg6+reg9; reg14=reg7+reg14; reg12=reg2+reg12; reg38=reg38+reg93;
   reg77=reg87+reg77; reg2=0.48803387171258487271*e.pos(6)[2]; reg45=reg82+reg45; reg37=reg37+reg96; reg7=reg11+reg36;
   reg69=reg69-reg85; double reg97=reg78-reg95; reg28=reg28-reg88; reg94=reg68+reg94; reg92=reg67+reg92;
   reg31=reg31-reg90; reg74=reg74-reg86; reg67=reg18+reg53; reg13=reg13-reg85; reg20=reg20+reg91;
   reg21=reg21-reg86; reg68=0.66666666666666670528*e.pos(7)[2]; reg24=reg24-reg87; reg22=reg22+reg89; reg12=reg17+reg12;
   reg28=reg28-reg86; reg59=reg59+reg7; reg97=reg82+reg97; reg31=reg31-reg85; reg77=reg77-reg68;
   reg75=reg22+reg75; reg22=0.66666666666666670528*e.pos(4)[0]; reg38=reg42+reg38; reg62=reg20+reg62; reg20=0.66666666666666670528*e.pos(4)[1];
   reg14=reg15+reg14; reg92=reg90+reg92; reg25=reg25+reg67; reg45=reg45-reg2; reg40=reg5+reg40;
   reg5=pow(reg69,2); reg42=pow(reg13,2); reg82=pow(reg21,2); reg33=reg3+reg33; reg37=reg30+reg37;
   reg24=reg24-reg68; reg3=pow(reg74,2); reg94=reg88+reg94; reg26=reg10+reg26; reg10=1.8213672050459180516*e.pos(5)[1];
   reg30=reg58+reg12; reg37=reg37-reg20; reg80=reg25+reg80; reg85=reg92-reg85; reg25=0.48803387171258487271*e.pos(5)[0];
   reg38=reg38-reg22; reg86=reg94-reg86; reg92=0.48803387171258487271*e.pos(5)[1]; reg94=0.66666666666666670528*e.pos(4)[2]; double reg98=reg50+reg14;
   reg45=reg45-reg68; reg40=reg51+reg40; reg82=reg42+reg82; reg33=reg44+reg33; reg42=pow(reg24,2);
   reg59=reg57+reg59; reg3=reg5+reg3; reg26=reg35+reg26; reg5=pow(reg28,2); reg52=reg19+reg52;
   reg19=pow(reg31,2); reg57=pow(reg77,2); reg75=reg75-reg22; double reg99=1.8213672050459180516*e.pos(5)[0]; reg97=reg97+reg2;
   reg62=reg62-reg20; double reg100=reg22+reg30; reg80=reg80-reg94; double reg101=0.48803387171258487271*e.pos(5)[2]; double reg102=pow(reg86,2);
   double reg103=0.66666666666666670528*e.pos(6)[1]; reg62=reg62-reg10; double reg104=pow(reg85,2); double reg105=0.66666666666666670528*e.pos(6)[0]; reg75=reg75-reg99;
   reg42=reg82+reg42; reg5=reg19+reg5; reg19=pow(reg45,2); reg68=reg97-reg68; reg82=reg20+reg98;
   reg38=reg38-reg25; reg57=reg3+reg57; reg52=reg65+reg52; reg37=reg37-reg92; reg3=reg43+reg33;
   reg97=reg71+reg26; reg59=reg59-reg94; double reg106=reg41+reg40; double reg107=1.8213672050459180516*e.pos(5)[2]; double reg108=pow(reg68,2);
   reg59=reg59-reg101; reg38=reg105+reg38; double reg109=1.8213672050459180516*e.pos(7)[0]; double reg110=reg55+reg52; double reg111=reg10-reg82;
   reg62=reg62+reg103; double reg112=0.48803387171258487271*e.pos(7)[1]; reg80=reg80-reg107; reg19=reg5+reg19; reg5=0.48803387171258487271*e.pos(7)[0];
   reg42=pow(reg42,0.5); reg22=reg22+reg106; reg102=reg104+reg102; reg75=reg75+reg105; reg57=pow(reg57,0.5);
   reg37=reg103+reg37; reg104=1.8213672050459180516*e.pos(7)[1]; reg20=reg20+reg3; double reg113=reg94+reg97; double reg114=reg99-reg100;
   double reg115=0.66666666666666670528*e.pos(6)[2]; double reg116=reg69/reg57; reg111=reg103+reg111; reg94=reg94+reg110; reg19=pow(reg19,0.5);
   double reg117=reg21/reg42; double reg118=reg13/reg42; double reg119=reg92-reg20; double reg120=reg25-reg22; double reg121=reg107-reg113;
   reg80=reg80+reg115; double reg122=0.48803387171258487271*e.pos(7)[2]; reg37=reg37-reg104; reg62=reg62-reg112; reg114=reg105+reg114;
   reg38=reg38-reg109; double reg123=1.8213672050459180516*e.pos(7)[2]; reg59=reg115+reg59; reg75=reg75-reg5; reg108=reg102+reg108;
   reg102=reg74/reg57; reg59=reg59-reg123; reg114=reg5+reg114; reg120=reg105+reg120; reg105=reg102*reg37;
   reg119=reg103+reg119; reg121=reg115+reg121; reg111=reg112+reg111; reg80=reg80-reg122; reg103=reg117*reg62;
   double reg124=reg116*reg38; double reg125=reg31/reg19; reg42=reg24/reg42; double reg126=reg28/reg19; double reg127=reg101-reg94;
   reg57=reg77/reg57; double reg128=reg118*reg75; reg108=pow(reg108,0.5); double reg129=reg126*reg111; reg19=reg45/reg19;
   double reg130=reg125*reg114; reg121=reg122+reg121; double reg131=reg42*reg80; reg105=reg124+reg105; reg124=reg57*reg59;
   reg127=reg115+reg127; reg115=reg85/reg108; reg103=reg128+reg103; reg128=reg86/reg108; reg119=reg104+reg119;
   reg120=reg109+reg120; reg131=reg103+reg131; reg124=reg105+reg124; reg103=reg128*reg119; reg105=reg115*reg120;
   reg108=reg68/reg108; double reg132=reg19*reg121; reg127=reg123+reg127; reg129=reg130+reg129; reg130=reg118*reg131;
   double reg133=reg117*reg131; reg103=reg105+reg103; reg105=reg116*reg124; double reg134=reg102*reg124; double reg135=reg108*reg127;
   reg132=reg129+reg132; reg129=reg125*reg132; double reg136=reg42*reg131; reg133=reg62-reg133; reg130=reg75-reg130;
   reg105=reg38-reg105; reg134=reg37-reg134; double reg137=reg57*reg124; double reg138=reg126*reg132; reg135=reg103+reg135;
   reg138=reg111-reg138; reg103=reg19*reg132; double reg139=reg128*reg135; double reg140=reg115*reg135; reg137=reg59-reg137;
   double reg141=pow(reg134,2); double reg142=pow(reg105,2); reg136=reg80-reg136; double reg143=pow(reg133,2); double reg144=pow(reg130,2);
   reg129=reg114-reg129; reg103=reg121-reg103; double reg145=pow(reg136,2); reg143=reg144+reg143; reg144=pow(reg137,2);
   double reg146=pow(reg138,2); double reg147=reg108*reg135; reg139=reg119-reg139; double reg148=pow(reg129,2); reg140=reg120-reg140;
   reg141=reg142+reg141; reg145=reg143+reg145; reg146=reg148+reg146; reg142=pow(reg140,2); reg143=pow(reg139,2);
   reg147=reg127-reg147; reg144=reg141+reg144; reg141=pow(reg103,2); reg145=pow(reg145,0.5); reg144=pow(reg144,0.5);
   reg148=pow(reg147,2); reg141=reg146+reg141; reg143=reg142+reg143; reg134=reg134/reg144; reg105=reg105/reg144;
   reg141=pow(reg141,0.5); reg133=reg133/reg145; reg130=reg130/reg145; reg148=reg143+reg148; reg142=reg13*reg130;
   reg143=reg21*reg133; reg116=reg69*reg116; reg102=reg74*reg102; reg148=pow(reg148,0.5); reg117=reg21*reg117;
   reg118=reg13*reg118; reg145=reg136/reg145; reg133=reg62*reg133; reg130=reg75*reg130; reg38=reg38*reg105;
   reg138=reg138/reg141; reg129=reg129/reg141; reg74=reg74*reg134; reg105=reg69*reg105; reg144=reg137/reg144;
   reg134=reg37*reg134; reg42=reg24*reg42; reg117=reg118+reg117; reg143=reg142+reg143; reg24=reg24*reg145;
   reg141=reg103/reg141; reg145=reg80*reg145; reg111=reg111*reg138; reg114=reg114*reg129; reg134=reg38+reg134;
   reg129=reg31*reg129; reg133=reg130+reg133; reg140=reg140/reg148; reg125=reg31*reg125; reg139=reg139/reg148;
   reg138=reg28*reg138; reg102=reg116+reg102; reg57=reg77*reg57; reg126=reg28*reg126; reg77=reg77*reg144;
   reg74=reg105+reg74; reg144=reg59*reg144; reg42=reg117+reg42; reg138=reg129+reg138; reg13=reg45*reg141;
   reg141=reg121*reg141; reg111=reg114+reg111; reg126=reg125+reg126; reg19=reg45*reg19; reg120=reg120*reg140;
   reg119=reg119*reg139; reg148=reg147/reg148; reg140=reg85*reg140; reg139=reg86*reg139; reg77=reg74+reg77;
   reg115=reg85*reg115; reg128=reg86*reg128; reg145=reg133+reg145; reg57=reg102+reg57; reg144=reg134+reg144;
   reg24=reg143+reg24; reg77=reg124*reg77; reg19=reg126+reg19; reg119=reg120+reg119; reg128=reg115+reg128;
   reg108=reg68*reg108; reg127=reg127*reg148; reg24=reg131*reg24; reg144=reg57*reg144; reg141=reg111+reg141;
   reg139=reg140+reg139; reg145=reg42*reg145; reg148=reg68*reg148; reg13=reg138+reg13; reg127=reg119+reg127;
   reg148=reg139+reg148; reg24=reg145-reg24; reg108=reg128+reg108; reg77=reg144-reg77; reg141=reg19*reg141;
   reg13=reg132*reg13; reg19=0.024056261216234395431*reg77; reg21=0.024056261216234395431*reg24; reg28=0.035220810900864524453*reg24; reg31=0.035220810900864524453*reg77;
   reg37=0.024056261216234409915*reg77; reg13=reg141-reg13; reg38=0.04166666666666666908*reg24; reg42=0.13144585576580215187*reg24; reg148=reg135*reg148;
   reg45=0.13144585576580215187*reg77; reg57=0.04166666666666666908*reg77; reg127=reg108*reg127; reg59=0.024056261216234409915*reg24; reg62=0.13144585576580215187*reg13;
   reg68=reg42+reg45; reg42=reg31+reg42; reg31=reg28+reg31; reg69=0.035220810900864524453*reg13; reg74=0.024056261216234409915*reg13;
   reg19=reg19-reg38; reg75=0.04166666666666666908*reg13; reg59=reg59+reg57; reg148=reg127-reg148; reg80=0.024056261216234395431*reg13;
   reg37=reg38+reg37; reg45=reg28+reg45; reg57=reg21-reg57; reg59=reg59+reg75; reg21=0.024056261216234409915*reg148;
   reg28=0.035220810900864524453*reg148; reg42=reg62+reg42; reg75=reg57-reg75; reg38=0.024056261216234395431*reg148; reg68=reg68+reg69;
   reg57=0.13144585576580215187*reg148; reg62=reg31+reg62; reg74=reg19-reg74; reg37=reg80-reg37; reg19=0.04166666666666666908*reg148;
   reg45=reg69+reg45; Ne(0,21)+=reg57+reg45; Ne(1,22)+=reg57+reg45; Ne(2,23)+=reg57+reg45; Ne(0,18)+=reg28+reg68;
   Ne(1,19)+=reg28+reg68; Ne(2,20)+=reg28+reg68; Ne(0,15)+=reg42+reg28; Ne(1,16)+=reg42+reg28; Ne(2,17)+=reg42+reg28;
   Ne(0,12)+=reg62+reg57; Ne(1,13)+=reg62+reg57; Ne(2,14)+=reg62+reg57; Ne(0,9)+=reg74-reg19; Ne(1,10)+=reg74-reg19;
   Ne(2,11)+=reg74-reg19; Ne(0,6)+=reg75-reg21; Ne(1,7)+=reg75-reg21; Ne(2,8)+=reg75-reg21; Ne(0,3)+=reg37-reg19;
   Ne(1,4)+=reg37-reg19; Ne(2,5)+=reg37-reg19; Ne(0,0)+=reg38-reg59; Ne(1,1)+=reg38-reg59; Ne(2,2)+=reg38-reg59;

}

};
} // namespace LMT

