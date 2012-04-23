
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
   
   double reg0=0.78379396366385999995*e.pos(1)[0]; double reg1=0.56758792732771999991*e.pos(0)[0]; double reg2=0.56758792732771999991*e.pos(1)[0]; double reg3=0.78379396366385999995*e.pos(0)[0]; double reg4=0.78379396366385999995*e.pos(1)[1];
   double reg5=0.56758792732771999991*e.pos(0)[1]; double reg6=0.78379396366385999995*e.pos(0)[1]; double reg7=0.56758792732771999991*e.pos(1)[1]; double reg8=1.3513818909915799999*e.pos(3)[1]; double reg9=reg4+reg5;
   double reg10=reg2+reg3; double reg11=0.78379396366385999995*e.pos(1)[2]; double reg12=1.3513818909915799999*e.pos(3)[0]; double reg13=reg0+reg1; double reg14=0.56758792732771999991*e.pos(0)[2];
   double reg15=reg7+reg6; double reg16=0.56758792732771999991*e.pos(1)[2]; double reg17=0.78379396366385999995*e.pos(0)[2]; double reg18=1.3513818909915799999*e.pos(3)[2]; double reg19=reg11+reg14;
   double reg20=0.63369514596091600003*e.pos(1)[0]; double reg21=reg12-reg10; double reg22=2.2673902919218320001*e.pos(0)[1]; double reg23=0.63369514596091600003*e.pos(1)[1]; double reg24=1.78379396366386*e.pos(4)[1];
   reg9=reg9-reg8; double reg25=2.2673902919218320001*e.pos(0)[0]; double reg26=reg8-reg15; double reg27=1.78379396366386*e.pos(4)[0]; reg13=reg13-reg12;
   double reg28=reg16+reg17; double reg29=2.9010854378827480001*e.pos(3)[1]; reg21=reg27+reg21; double reg30=2.2673902919218320001*e.pos(0)[2]; double reg31=0.63369514596091600003*e.pos(1)[2];
   double reg32=reg23+reg22; double reg33=0.63369514596091600003*e.pos(0)[1]; double reg34=2.2673902919218320001*e.pos(1)[1]; double reg35=0.63369514596091600003*e.pos(0)[0]; reg4=reg4-reg6;
   reg0=reg0-reg3; double reg36=0.43241207267228000009*e.pos(4)[1]; double reg37=2.2673902919218320001*e.pos(1)[0]; double reg38=reg18-reg28; double reg39=0.43241207267228000009*e.pos(4)[0];
   double reg40=reg20+reg25; reg13=reg13+reg27; double reg41=1.78379396366386*e.pos(5)[0]; reg9=reg9+reg24; double reg42=1.78379396366386*e.pos(5)[1];
   double reg43=2.9010854378827480001*e.pos(3)[0]; reg26=reg24+reg26; reg19=reg19-reg18; double reg44=1.78379396366386*e.pos(4)[2]; double reg45=reg29-reg32;
   double reg46=reg43-reg40; double reg47=0.366304854039084*e.pos(4)[1]; reg4=reg36+reg4; double reg48=0.366304854039084*e.pos(4)[0]; reg11=reg11-reg17;
   double reg49=0.43241207267228000009*e.pos(5)[1]; reg38=reg44+reg38; reg13=reg13-reg41; reg26=reg26-reg42; reg42=reg9-reg42;
   reg19=reg19+reg44; reg9=1.78379396366386*e.pos(5)[2]; double reg50=reg31+reg30; double reg51=2.9010854378827480001*e.pos(3)[2]; double reg52=0.43241207267228000009*e.pos(5)[0];
   double reg53=0.43241207267228000009*e.pos(4)[2]; reg0=reg39+reg0; reg37=reg37+reg35; reg34=reg34+reg33; double reg54=2.2673902919218320001*e.pos(1)[2];
   double reg55=0.63369514596091600003*e.pos(0)[2]; reg41=reg21-reg41; reg11=reg53+reg11; reg21=pow(reg26,2); reg49=reg4-reg49;
   reg52=reg0-reg52; reg0=pow(reg41,2); reg38=reg38-reg9; reg4=0.43241207267228000009*e.pos(5)[2]; reg34=reg34-reg29;
   double reg56=0.78379396366385999995*e.pos(2)[0]; double reg57=pow(reg42,2); reg37=reg37-reg43; double reg58=pow(reg13,2); double reg59=2.710505431213761085e-20*e.pos(3)[0];
   double reg60=0.78379396366385999995*e.pos(2)[1]; double reg61=0.366304854039084*e.pos(4)[2]; double reg62=reg51-reg50; reg9=reg19-reg9; reg19=2.710505431213761085e-20*e.pos(3)[1];
   double reg63=0.366304854039084*e.pos(5)[1]; reg45=reg45+reg47; double reg64=reg33-reg23; double reg65=0.366304854039084*e.pos(5)[0]; reg54=reg54+reg55;
   reg46=reg46+reg48; double reg66=reg35-reg20; double reg67=3.2673902919218320001*e.pos(4)[1]; reg66=reg66-reg59; reg64=reg64-reg19;
   double reg68=reg55-reg31; double reg69=2.710505431213761085e-20*e.pos(3)[2]; double reg70=3.2673902919218320001*e.pos(4)[0]; reg54=reg54-reg51; reg34=reg47+reg34;
   reg37=reg48+reg37; double reg71=0.366304854039084*e.pos(5)[2]; reg62=reg62+reg61; reg45=reg45-reg63; double reg72=pow(reg52,2);
   reg46=reg46-reg65; double reg73=pow(reg49,2); double reg74=0.56758792732771999991*e.pos(2)[1]; reg4=reg11-reg4; reg11=0.56758792732771999991*e.pos(2)[0];
   double reg75=reg1+reg56; double reg76=1.78379396366386*e.pos(3)[0]; reg21=reg0+reg21; reg0=pow(reg9,2); reg57=reg58+reg57;
   reg58=reg5+reg60; double reg77=1.78379396366386*e.pos(3)[1]; double reg78=0.78379396366385999995*e.pos(2)[2]; double reg79=pow(reg38,2); double reg80=3.2673902919218320001*e.pos(5)[0];
   reg64=reg67+reg64; reg66=reg70+reg66; double reg81=3.2673902919218320001*e.pos(5)[1]; double reg82=3.2673902919218320001*e.pos(4)[2]; reg0=reg57+reg0;
   reg68=reg68-reg69; reg73=reg72+reg73; reg54=reg61+reg54; reg57=pow(reg4,2); reg63=reg34-reg63;
   reg34=reg3+reg11; reg65=reg37-reg65; reg37=reg6+reg74; reg75=reg75-reg76; reg72=0.56758792732771999991*e.pos(2)[2];
   double reg83=0.63369514596091600003*e.pos(2)[1]; double reg84=0.63369514596091600003*e.pos(2)[0]; reg62=reg62-reg71; reg58=reg58-reg77; double reg85=pow(reg45,2);
   double reg86=reg14+reg78; double reg87=1.78379396366386*e.pos(3)[2]; double reg88=pow(reg46,2); reg79=reg21+reg79; reg21=0.43241207267228000009*e.pos(3)[1];
   reg6=reg60-reg6; reg3=reg56-reg3; reg56=0.43241207267228000009*e.pos(3)[0]; reg71=reg54-reg71; reg57=reg73+reg57;
   reg54=pow(reg63,2); reg21=reg6-reg21; reg76=reg76+reg34; reg6=pow(reg65,2); reg77=reg77+reg37;
   reg56=reg3-reg56; reg3=reg17+reg72; reg75=reg27+reg75; reg60=1.3513818909915799999*e.pos(5)[0]; reg73=0.63369514596091600003*e.pos(2)[2];
   double reg89=0.366304854039084*e.pos(3)[1]; reg22=reg22+reg83; double reg90=0.366304854039084*e.pos(3)[0]; reg25=reg25+reg84; double reg91=pow(reg62,2);
   reg58=reg24+reg58; double reg92=1.3513818909915799999*e.pos(5)[1]; reg85=reg88+reg85; reg79=pow(reg79,0.5); reg86=reg86-reg87;
   reg17=reg78-reg17; reg78=0.43241207267228000009*e.pos(3)[2]; reg80=reg66-reg80; reg81=reg64-reg81; reg0=pow(reg0,0.5);
   reg68=reg82+reg68; reg64=3.2673902919218320001*e.pos(5)[2]; reg24=reg24-reg77; reg39=reg56+reg39; reg87=reg87+reg3;
   reg56=1.3513818909915799999*e.pos(5)[2]; reg66=0.366304854039084*e.pos(3)[2]; reg30=reg30+reg73; reg78=reg17-reg78; reg17=reg22+reg89;
   reg88=2.2673902919218320001*e.pos(2)[0]; double reg93=reg25+reg90; reg86=reg44+reg86; double reg94=pow(reg80,2); double reg95=2.2673902919218320001*e.pos(2)[1];
   reg91=reg85+reg91; reg75=reg75-reg60; reg64=reg68-reg64; reg58=reg58-reg92; reg68=reg41/reg79;
   reg85=pow(reg81,2); double reg96=reg26/reg79; reg57=pow(reg57,0.5); double reg97=reg33-reg83; reg36=reg21+reg36;
   reg21=3.2673902919218320001*e.pos(3)[1]; double reg98=pow(reg71,2); reg27=reg27-reg76; double reg99=reg13/reg0; double reg100=reg42/reg0;
   double reg101=reg35-reg84; double reg102=3.2673902919218320001*e.pos(3)[0]; reg54=reg6+reg54; reg91=pow(reg91,0.5); reg6=reg49/reg57;
   double reg103=reg52/reg57; double reg104=pow(reg64,2); reg95=reg33+reg95; double reg105=reg48-reg93; double reg106=reg99*reg75;
   reg79=reg38/reg79; double reg107=2.2673902919218320001*e.pos(2)[2]; double reg108=reg96*reg36; reg98=reg54+reg98; reg54=reg100*reg58;
   reg88=reg35+reg88; reg85=reg94+reg85; reg86=reg86-reg56; reg102=reg101-reg102; reg94=reg68*reg39;
   reg21=reg97-reg21; reg24=reg92+reg24; reg27=reg60+reg27; reg97=reg55-reg73; reg101=3.2673902919218320001*e.pos(3)[2];
   reg44=reg44-reg87; double reg109=reg30+reg66; reg0=reg9/reg0; double reg110=2.9010854378827480001*e.pos(5)[0]; reg53=reg78+reg53;
   reg78=2.9010854378827480001*e.pos(5)[1]; double reg111=reg47-reg17; double reg112=5.42101086242752217e-20*e.pos(5)[0]; reg67=reg21+reg67; reg105=reg105+reg110;
   reg104=reg85+reg104; reg98=pow(reg98,0.5); reg111=reg111+reg78; reg21=5.42101086242752217e-20*e.pos(5)[1]; reg107=reg55+reg107;
   reg54=reg106+reg54; reg90=reg88-reg90; reg85=reg45/reg91; reg88=reg46/reg91; reg101=reg97-reg101;
   reg97=2.9010854378827480001*e.pos(5)[2]; reg70=reg102+reg70; reg102=reg61-reg109; reg89=reg95-reg89; reg95=reg0*reg86;
   reg108=reg94+reg108; reg57=reg4/reg57; reg94=reg103*reg27; reg106=reg6*reg24; reg44=reg56+reg44;
   double reg113=reg79*reg53; reg66=reg107-reg66; reg91=reg62/reg91; reg89=reg47+reg89; reg90=reg48+reg90;
   reg47=reg88*reg105; reg48=reg85*reg111; reg102=reg102+reg97; reg104=pow(reg104,0.5); reg113=reg108+reg113;
   reg107=reg57*reg44; reg106=reg94+reg106; reg94=reg65/reg98; reg108=reg63/reg98; double reg114=5.42101086242752217e-20*e.pos(5)[2];
   reg82=reg101+reg82; reg67=reg67-reg21; reg70=reg70-reg112; reg95=reg54+reg95; reg54=reg68*reg113;
   reg101=reg80/reg104; double reg115=reg96*reg113; double reg116=reg81/reg104; reg48=reg47+reg48; reg47=reg100*reg95;
   reg90=reg90-reg110; double reg117=reg91*reg102; reg98=reg71/reg98; reg89=reg89-reg78; reg82=reg82-reg114;
   reg107=reg106+reg107; reg106=reg108*reg67; double reg118=reg99*reg95; double reg119=reg94*reg70; reg66=reg61+reg66;
   reg61=reg103*reg107; reg117=reg48+reg117; reg104=reg64/reg104; reg54=reg39-reg54; reg115=reg36-reg115;
   reg48=reg0*reg95; double reg120=reg79*reg113; reg47=reg58-reg47; double reg121=reg98*reg82; reg118=reg75-reg118;
   reg106=reg119+reg106; reg66=reg66-reg97; reg119=reg6*reg107; double reg122=reg116*reg89; double reg123=reg101*reg90;
   double reg124=reg104*reg66; reg122=reg123+reg122; reg123=pow(reg118,2); reg121=reg106+reg121; reg106=pow(reg47,2);
   reg120=reg53-reg120; double reg125=pow(reg115,2); reg48=reg86-reg48; double reg126=pow(reg54,2); double reg127=reg57*reg107;
   reg61=reg27-reg61; reg119=reg24-reg119; double reg128=reg88*reg117; double reg129=reg85*reg117; double reg130=pow(reg61,2);
   reg124=reg122+reg124; reg128=reg105-reg128; reg122=reg94*reg121; double reg131=reg108*reg121; reg126=reg125+reg126;
   reg106=reg123+reg106; reg129=reg111-reg129; reg123=pow(reg48,2); reg125=pow(reg120,2); double reg132=reg91*reg117;
   reg127=reg44-reg127; double reg133=pow(reg119,2); double reg134=pow(reg128,2); double reg135=reg98*reg121; reg131=reg67-reg131;
   double reg136=reg116*reg124; double reg137=pow(reg129,2); double reg138=pow(reg127,2); reg123=reg106+reg123; reg122=reg70-reg122;
   reg133=reg130+reg133; reg106=reg101*reg124; reg125=reg126+reg125; reg132=reg102-reg132; reg106=reg90-reg106;
   reg126=pow(reg122,2); reg136=reg89-reg136; reg130=pow(reg132,2); double reg139=pow(reg131,2); reg135=reg82-reg135;
   reg137=reg134+reg137; reg138=reg133+reg138; reg133=reg104*reg124; reg125=pow(reg125,0.5); reg123=pow(reg123,0.5);
   reg130=reg137+reg130; reg134=pow(reg136,2); reg133=reg66-reg133; reg54=reg54/reg125; reg115=reg115/reg125;
   reg137=pow(reg106,2); reg138=pow(reg138,0.5); reg47=reg47/reg123; reg118=reg118/reg123; reg139=reg126+reg139;
   reg126=pow(reg135,2); reg36=reg36*reg115; double reg140=pow(reg133,2); reg39=reg39*reg54; reg125=reg120/reg125;
   reg134=reg137+reg134; reg115=reg26*reg115; reg96=reg26*reg96; reg68=reg41*reg68; reg126=reg139+reg126;
   reg54=reg41*reg54; reg123=reg48/reg123; reg58=reg58*reg47; reg99=reg13*reg99; reg100=reg42*reg100;
   reg130=pow(reg130,0.5); reg61=reg61/reg138; reg75=reg75*reg118; reg47=reg42*reg47; reg119=reg119/reg138;
   reg118=reg13*reg118; reg58=reg75+reg58; reg103=reg52*reg103; reg6=reg49*reg6; reg126=pow(reg126,0.5);
   reg79=reg38*reg79; reg24=reg24*reg119; reg27=reg27*reg61; reg128=reg128/reg130; reg129=reg129/reg130;
   reg36=reg39+reg36; reg0=reg9*reg0; reg140=reg134+reg140; reg100=reg99+reg100; reg53=reg53*reg125;
   reg119=reg49*reg119; reg61=reg52*reg61; reg9=reg9*reg123; reg47=reg118+reg47; reg125=reg38*reg125;
   reg115=reg54+reg115; reg138=reg127/reg138; reg96=reg68+reg96; reg123=reg86*reg123; reg53=reg36+reg53;
   reg13=reg45*reg129; reg9=reg47+reg9; reg26=reg46*reg128; reg36=reg4*reg138; reg57=reg4*reg57;
   reg140=pow(reg140,0.5); reg119=reg61+reg119; reg130=reg132/reg130; reg129=reg111*reg129; reg123=reg58+reg123;
   reg128=reg105*reg128; reg88=reg46*reg88; reg125=reg115+reg125; reg138=reg44*reg138; reg24=reg27+reg24;
   reg79=reg96+reg79; reg131=reg131/reg126; reg0=reg100+reg0; reg122=reg122/reg126; reg85=reg45*reg85;
   reg6=reg103+reg6; reg13=reg26+reg13; reg4=reg62*reg130; reg130=reg102*reg130; reg129=reg128+reg129;
   reg91=reg62*reg91; reg85=reg88+reg85; reg138=reg24+reg138; reg36=reg119+reg36; reg126=reg135/reg126;
   reg67=reg67*reg131; reg70=reg70*reg122; reg94=reg65*reg94; reg122=reg65*reg122; reg131=reg63*reg131;
   reg108=reg63*reg108; reg123=reg0*reg123; reg136=reg136/reg140; reg125=reg113*reg125; reg106=reg106/reg140;
   reg53=reg79*reg53; reg57=reg6+reg57; reg9=reg95*reg9; reg125=reg53-reg125; reg138=reg57*reg138;
   reg98=reg71*reg98; reg91=reg85+reg91; reg108=reg94+reg108; reg9=reg123-reg9; reg36=reg107*reg36;
   reg67=reg70+reg67; reg4=reg13+reg4; reg0=reg81*reg136; reg6=reg80*reg106; reg130=reg129+reg130;
   reg106=reg90*reg106; reg136=reg89*reg136; reg101=reg80*reg101; reg71=reg71*reg126; reg116=reg81*reg116;
   reg131=reg122+reg131; reg140=reg133/reg140; reg126=reg82*reg126; reg66=reg66*reg140; reg136=reg106+reg136;
   reg13=0.009463616120767210603*reg125; reg24=0.088847818743090689935*reg9; reg26=0.088847818743090689935*reg125; reg27=0.005384432036113586778*reg9; reg140=reg64*reg140;
   reg38=0.021537728144454347112*reg125; reg39=0.021537728144454347112*reg9; reg0=reg6+reg0; reg104=reg64*reg104; reg116=reg101+reg116;
   reg6=0.009463616120767210603*reg9; reg4=reg117*reg4; reg130=reg91*reg130; reg36=reg138-reg36; reg71=reg131+reg71;
   reg41=0.005384432036113586778*reg125; reg126=reg67+reg126; reg98=reg108+reg98; reg26=reg39+reg26; reg13=reg27+reg13;
   reg27=reg41+reg27; reg42=0.009463616120767210603*reg36; reg44=0.021537728144454347112*reg36; reg24=reg38+reg24; reg45=0.005384432036113586778*reg36;
   reg38=reg39+reg38; reg41=reg6+reg41; reg6=0.088847818743090689935*reg36; reg104=reg116+reg104; reg140=reg0+reg140;
   reg4=reg130-reg4; reg71=reg121*reg71; reg126=reg98*reg126; reg66=reg136+reg66; reg0=0.016449618187943419918*reg4;
   reg6=reg38+reg6; reg71=reg126-reg71; reg41=reg41+reg45; reg66=reg104*reg66; reg42=reg27+reg42;
   reg24=reg24+reg44; reg27=0.028457289286966203713*reg4; reg38=0.0018441552587796664112*reg4; reg39=0.0041124045469858549794*reg4; reg13=reg45+reg13;
   reg140=reg124*reg140; reg26=reg44+reg26; reg26=reg0+reg26; reg38=reg24+reg38; reg24=0.0018441552587796664109*reg71;
   reg44=0.016449618187943419918*reg71; reg45=0.016449618187943419916*reg71; reg46=0.004112404546985854979*reg71; reg41=reg27-reg41; reg140=reg66-reg140;
   reg27=0.028457289286966203713*reg71; reg13=reg13+reg39; reg42=reg39+reg42; reg39=0.0041124045469858549794*reg71; reg0=reg6+reg0;
   reg6=0.028457289286966203713*reg140; reg39=reg42+reg39; reg13=reg27-reg13; reg24=reg26+reg24; reg26=0.0018441552587796664111*reg140;
   reg44=reg38+reg44; reg27=0.016449618187943419918*reg140; reg46=reg41-reg46; reg38=0.0041124045469858549794*reg140; reg45=reg0+reg45;
   Ne(0,15)+=reg27+reg24; Ne(1,16)+=reg27+reg24; Ne(2,17)+=reg27+reg24; Ne(0,6)+=reg6-reg39; Ne(1,7)+=reg6-reg39;
   Ne(2,8)+=reg6-reg39; Ne(0,0)+=reg46-reg38; Ne(1,1)+=reg46-reg38; Ne(2,2)+=reg46-reg38; Ne(0,12)+=reg44+reg27;
   Ne(1,13)+=reg44+reg27; Ne(2,14)+=reg44+reg27; Ne(0,9)+=reg45+reg26; Ne(1,10)+=reg45+reg26; Ne(2,11)+=reg45+reg26;
   Ne(0,3)+=reg13-reg38; Ne(1,4)+=reg13-reg38; Ne(2,5)+=reg13-reg38;

}
 /// pour les Quad_8
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<Quad_8,TNB,TN,TD,NET> e,TMAT &Ne) const {
   
   double reg0=0.36602540378443865451*e.pos(1)[1]; double reg1=0.12200846792814621817*e.pos(0)[1]; double reg2=0.12200846792814621817*e.pos(1)[0]; double reg3=0.36602540378443865451*e.pos(0)[0]; double reg4=0.36602540378443865451*e.pos(0)[1];
   double reg5=0.12200846792814621817*e.pos(1)[1]; double reg6=0.36602540378443865451*e.pos(1)[0]; double reg7=0.12200846792814621817*e.pos(0)[0]; double reg8=reg5+reg4; double reg9=1.3660254037844385386*e.pos(2)[1];
   double reg10=0.45534180126147951289*e.pos(0)[0]; double reg11=0.45534180126147951289*e.pos(0)[1]; double reg12=0.12200846792814621817*e.pos(0)[2]; double reg13=1.3660254037844385386*e.pos(1)[0]; double reg14=0.36602540378443865451*e.pos(1)[2];
   double reg15=reg2+reg3; double reg16=0.45534180126147951289*e.pos(2)[0]; double reg17=0.12200846792814621817*e.pos(1)[2]; double reg18=reg6+reg7; double reg19=1.3660254037844385386*e.pos(2)[0];
   double reg20=0.45534180126147951289*e.pos(2)[1]; double reg21=reg0+reg1; double reg22=0.36602540378443865451*e.pos(0)[2]; double reg23=1.3660254037844385386*e.pos(1)[1]; double reg24=0.45534180126147951289*e.pos(1)[0];
   double reg25=0.45534180126147951289*e.pos(1)[1]; double reg26=reg17+reg22; reg8=reg8+reg9; double reg27=1.3660254037844385386*e.pos(2)[2]; double reg28=0.45534180126147951289*e.pos(3)[1];
   double reg29=0.45534180126147951289*e.pos(2)[2]; double reg30=reg14+reg12; double reg31=1.3660254037844385386*e.pos(3)[0]; double reg32=1.3660254037844385386*e.pos(3)[1]; double reg33=reg21+reg20;
   double reg34=reg16+reg18; double reg35=reg11+reg23; double reg36=1.3660254037844385386*e.pos(1)[2]; double reg37=reg10+reg13; reg15=reg19+reg15;
   double reg38=1.3660254037844385386*e.pos(0)[1]; double reg39=0.45534180126147951289*e.pos(3)[0]; double reg40=0.45534180126147951289*e.pos(0)[2]; double reg41=0.12200846792814621817*e.pos(2)[1]; double reg42=0.12200846792814621817*e.pos(2)[0];
   double reg43=1.3660254037844385386*e.pos(0)[0]; double reg44=reg24+reg43; double reg45=reg25+reg38; double reg46=reg40+reg36; double reg47=0.45534180126147951289*e.pos(1)[2];
   double reg48=reg34+reg31; double reg49=0.36602540378443865451*e.pos(2)[0]; double reg50=0.36602540378443865451*e.pos(3)[1]; reg35=reg41+reg35; double reg51=0.36602540378443865451*e.pos(3)[0];
   reg37=reg42+reg37; double reg52=0.36602540378443865451*e.pos(2)[1]; double reg53=reg33+reg32; double reg54=1.3660254037844385386*e.pos(0)[2]; double reg55=0.12200846792814621817*e.pos(2)[2];
   double reg56=1.3660254037844385386*e.pos(3)[2]; double reg57=reg30+reg29; double reg58=0.48803387171258487271*e.pos(4)[0]; reg15=reg15+reg39; double reg59=0.45534180126147951289*e.pos(3)[2];
   reg26=reg26+reg27; double reg60=0.48803387171258487271*e.pos(4)[1]; reg8=reg8+reg28; reg35=reg35+reg50; double reg61=1.8213672050459180516*e.pos(4)[1];
   double reg62=1.8213672050459180516*e.pos(4)[0]; reg37=reg37+reg51; double reg63=0.66666666666666670528*e.pos(5)[0]; double reg64=reg45+reg52; double reg65=reg47+reg54;
   double reg66=0.36602540378443865451*e.pos(2)[2]; reg15=reg15-reg58; double reg67=0.48803387171258487271*e.pos(4)[2]; reg26=reg26+reg59; double reg68=reg58-reg48;
   double reg69=reg57+reg56; double reg70=reg60-reg53; double reg71=0.66666666666666670528*e.pos(5)[1]; reg8=reg8-reg60; double reg72=0.12200846792814621817*e.pos(3)[0];
   double reg73=0.36602540378443865451*e.pos(3)[2]; reg46=reg55+reg46; double reg74=0.12200846792814621817*e.pos(3)[1]; double reg75=reg49+reg44; reg37=reg37-reg62;
   reg46=reg46+reg73; double reg76=1.8213672050459180516*e.pos(4)[2]; double reg77=reg74+reg64; double reg78=reg65+reg66; reg15=reg63+reg15;
   double reg79=1.8213672050459180516*e.pos(6)[0]; double reg80=0.66666666666666670528*e.pos(5)[2]; double reg81=reg72+reg75; reg70=reg71+reg70; reg26=reg26-reg67;
   reg68=reg63+reg68; double reg82=reg67-reg69; double reg83=1.8213672050459180516*e.pos(6)[1]; reg8=reg8+reg71; double reg84=0.12200846792814621817*e.pos(3)[2];
   reg35=reg35-reg61; double reg85=reg3+reg24; double reg86=reg0+reg11; double reg87=1.8213672050459180516*e.pos(6)[2]; double reg88=reg6+reg10;
   reg15=reg15-reg79; double reg89=0.66666666666666670528*e.pos(7)[0]; reg26=reg26+reg80; reg68=reg79+reg68; reg82=reg80+reg82;
   reg35=reg71+reg35; double reg90=0.66666666666666670528*e.pos(7)[1]; reg8=reg8-reg83; double reg91=0.48803387171258487271*e.pos(6)[1]; reg70=reg83+reg70;
   double reg92=reg62-reg81; double reg93=reg4+reg25; reg37=reg63+reg37; double reg94=0.48803387171258487271*e.pos(6)[0]; double reg95=reg61-reg77;
   double reg96=reg84+reg78; reg46=reg46-reg76; reg13=reg7+reg13; reg23=reg1+reg23; reg68=reg68-reg89;
   reg1=0.48803387171258487271*e.pos(6)[2]; reg8=reg8-reg90; reg46=reg80+reg46; reg35=reg35-reg91; reg92=reg63+reg92;
   reg37=reg37-reg94; reg70=reg70-reg90; reg15=reg15-reg89; reg9=reg9+reg93; reg7=reg14+reg40;
   reg63=reg22+reg47; reg41=reg41+reg86; double reg97=0.66666666666666670528*e.pos(7)[2]; reg42=reg42+reg88; double reg98=reg76-reg96;
   reg19=reg19+reg85; reg82=reg87+reg82; reg26=reg26-reg87; reg95=reg71+reg95; reg92=reg94+reg92;
   reg13=reg16+reg13; reg41=reg32+reg41; reg32=pow(reg15,2); reg55=reg55+reg7; reg23=reg20+reg23;
   reg98=reg80+reg98; reg35=reg35-reg90; reg42=reg31+reg42; reg37=reg37-reg89; reg95=reg91+reg95;
   reg43=reg2+reg43; reg82=reg82-reg97; reg36=reg12+reg36; reg46=reg46-reg1; reg38=reg5+reg38;
   reg27=reg27+reg63; reg2=0.66666666666666670528*e.pos(4)[1]; reg74=reg9+reg74; reg5=0.66666666666666670528*e.pos(4)[0]; reg72=reg19+reg72;
   reg9=pow(reg68,2); reg12=pow(reg8,2); reg26=reg26-reg97; reg19=pow(reg70,2); reg41=reg41-reg2;
   reg31=0.48803387171258487271*e.pos(5)[1]; reg54=reg17+reg54; reg55=reg56+reg55; reg98=reg1+reg98; reg90=reg95-reg90;
   reg17=0.66666666666666670528*e.pos(4)[2]; reg56=pow(reg37,2); reg84=reg27+reg84; reg27=1.8213672050459180516*e.pos(5)[1]; reg74=reg74-reg2;
   reg71=pow(reg35,2); reg80=pow(reg26,2); reg95=1.8213672050459180516*e.pos(5)[0]; reg72=reg72-reg5; reg46=reg46-reg97;
   reg36=reg29+reg36; reg89=reg92-reg89; reg92=reg51+reg13; double reg99=reg50+reg23; reg43=reg49+reg43;
   double reg100=pow(reg82,2); reg38=reg52+reg38; reg42=reg42-reg5; double reg101=0.48803387171258487271*e.pos(5)[0]; reg12=reg32+reg12;
   reg19=reg9+reg19; reg42=reg42-reg101; reg97=reg98-reg97; reg9=1.8213672050459180516*e.pos(5)[2]; reg84=reg84-reg17;
   reg32=reg28+reg38; reg98=0.66666666666666670528*e.pos(6)[1]; reg74=reg74-reg27; double reg102=pow(reg90,2); reg71=reg56+reg71;
   reg56=0.66666666666666670528*e.pos(6)[0]; reg72=reg72-reg95; double reg103=reg39+reg43; double reg104=pow(reg46,2); double reg105=0.48803387171258487271*e.pos(5)[2];
   reg55=reg55-reg17; double reg106=pow(reg89,2); reg54=reg66+reg54; reg80=reg12+reg80; reg100=reg19+reg100;
   reg12=reg5+reg92; reg19=reg2+reg99; double reg107=reg73+reg36; reg41=reg41-reg31; reg102=reg106+reg102;
   reg106=0.66666666666666670528*e.pos(6)[2]; reg84=reg84-reg9; reg80=pow(reg80,0.5); double reg108=reg17+reg107; reg5=reg5+reg103;
   double reg109=reg27-reg19; double reg110=0.48803387171258487271*e.pos(7)[1]; reg74=reg74+reg98; double reg111=reg95-reg12; double reg112=0.48803387171258487271*e.pos(7)[0];
   reg72=reg72+reg56; reg100=pow(reg100,0.5); reg104=reg71+reg104; reg2=reg2+reg32; reg41=reg98+reg41;
   reg71=1.8213672050459180516*e.pos(7)[1]; double reg113=1.8213672050459180516*e.pos(7)[0]; double reg114=reg59+reg54; reg42=reg56+reg42; reg55=reg55-reg105;
   double reg115=pow(reg97,2); double reg116=reg9-reg108; double reg117=reg68/reg100; reg104=pow(reg104,0.5); reg55=reg106+reg55;
   double reg118=1.8213672050459180516*e.pos(7)[2]; double reg119=reg70/reg100; reg111=reg56+reg111; reg72=reg72-reg112; reg42=reg42-reg113;
   double reg120=reg101-reg5; double reg121=reg15/reg80; reg74=reg74-reg110; reg41=reg41-reg71; reg109=reg98+reg109;
   double reg122=reg8/reg80; reg115=reg102+reg115; reg102=reg31-reg2; reg17=reg17+reg114; reg84=reg84+reg106;
   double reg123=0.48803387171258487271*e.pos(7)[2]; reg116=reg106+reg116; reg80=reg26/reg80; double reg124=reg35/reg104; double reg125=reg37/reg104;
   reg102=reg98+reg102; reg55=reg55-reg118; reg98=reg119*reg41; reg84=reg84-reg123; reg111=reg112+reg111;
   double reg126=reg121*reg72; reg120=reg56+reg120; reg100=reg82/reg100; reg56=reg117*reg42; reg115=pow(reg115,0.5);
   double reg127=reg105-reg17; double reg128=reg122*reg74; reg109=reg110+reg109; reg102=reg71+reg102; reg128=reg126+reg128;
   reg104=reg46/reg104; reg127=reg106+reg127; reg116=reg123+reg116; reg106=reg89/reg115; reg98=reg56+reg98;
   reg56=reg124*reg109; reg126=reg80*reg84; reg120=reg113+reg120; double reg129=reg100*reg55; double reg130=reg90/reg115;
   double reg131=reg125*reg111; reg56=reg131+reg56; reg127=reg118+reg127; reg131=reg104*reg116; reg129=reg98+reg129;
   reg126=reg128+reg126; reg98=reg130*reg102; reg115=reg97/reg115; reg128=reg106*reg120; double reg132=reg121*reg126;
   double reg133=reg122*reg126; double reg134=reg115*reg127; double reg135=reg119*reg129; reg131=reg56+reg131; reg56=reg117*reg129;
   reg98=reg128+reg98; reg128=reg124*reg131; double reg136=reg125*reg131; double reg137=reg80*reg126; reg132=reg72-reg132;
   reg134=reg98+reg134; reg98=reg100*reg129; reg133=reg74-reg133; reg56=reg42-reg56; reg135=reg41-reg135;
   double reg138=pow(reg132,2); reg98=reg55-reg98; double reg139=reg104*reg131; double reg140=pow(reg135,2); reg128=reg109-reg128;
   reg136=reg111-reg136; double reg141=pow(reg56,2); double reg142=pow(reg133,2); double reg143=reg106*reg134; double reg144=reg130*reg134;
   reg137=reg84-reg137; reg144=reg102-reg144; reg143=reg120-reg143; double reg145=reg115*reg134; double reg146=pow(reg136,2);
   double reg147=pow(reg128,2); reg139=reg116-reg139; double reg148=pow(reg137,2); reg142=reg138+reg142; reg138=pow(reg98,2);
   reg140=reg141+reg140; reg147=reg146+reg147; reg141=pow(reg139,2); reg138=reg140+reg138; reg148=reg142+reg148;
   reg145=reg127-reg145; reg140=pow(reg143,2); reg142=pow(reg144,2); reg148=pow(reg148,0.5); reg142=reg140+reg142;
   reg138=pow(reg138,0.5); reg141=reg147+reg141; reg140=pow(reg145,2); reg132=reg132/reg148; reg133=reg133/reg148;
   reg56=reg56/reg138; reg141=pow(reg141,0.5); reg135=reg135/reg138; reg140=reg142+reg140; reg42=reg42*reg56;
   reg41=reg41*reg135; reg138=reg98/reg138; reg119=reg70*reg119; reg117=reg68*reg117; reg98=reg8*reg133;
   reg142=reg15*reg132; reg148=reg137/reg148; reg133=reg74*reg133; reg132=reg72*reg132; reg122=reg8*reg122;
   reg121=reg15*reg121; reg140=pow(reg140,0.5); reg136=reg136/reg141; reg128=reg128/reg141; reg135=reg70*reg135;
   reg56=reg68*reg56; reg84=reg84*reg148; reg133=reg132+reg133; reg98=reg142+reg98; reg148=reg26*reg148;
   reg8=reg35*reg128; reg15=reg37*reg136; reg80=reg26*reg80; reg122=reg121+reg122; reg124=reg35*reg124;
   reg119=reg117+reg119; reg141=reg139/reg141; reg100=reg82*reg100; reg128=reg109*reg128; reg55=reg55*reg138;
   reg41=reg42+reg41; reg135=reg56+reg135; reg138=reg82*reg138; reg144=reg144/reg140; reg125=reg37*reg125;
   reg143=reg143/reg140; reg136=reg111*reg136; reg120=reg120*reg143; reg116=reg116*reg141; reg102=reg102*reg144;
   reg140=reg145/reg140; reg80=reg122+reg80; reg143=reg89*reg143; reg128=reg136+reg128; reg8=reg15+reg8;
   reg141=reg46*reg141; reg144=reg90*reg144; reg55=reg41+reg55; reg138=reg135+reg138; reg106=reg89*reg106;
   reg100=reg119+reg100; reg130=reg90*reg130; reg124=reg125+reg124; reg148=reg98+reg148; reg104=reg46*reg104;
   reg84=reg133+reg84; reg102=reg120+reg102; reg138=reg129*reg138; reg84=reg80*reg84; reg127=reg127*reg140;
   reg55=reg100*reg55; reg130=reg106+reg130; reg115=reg97*reg115; reg144=reg143+reg144; reg140=reg97*reg140;
   reg116=reg128+reg116; reg148=reg126*reg148; reg141=reg8+reg141; reg104=reg124+reg104; reg140=reg144+reg140;
   reg127=reg102+reg127; reg138=reg55-reg138; reg148=reg84-reg148; reg115=reg130+reg115; reg141=reg131*reg141;
   reg116=reg104*reg116; reg8=0.024056261216234395431*reg148; reg15=0.13144585576580215187*reg138; reg26=0.024056261216234395431*reg138; reg35=0.024056261216234409915*reg138;
   reg37=0.13144585576580215187*reg148; reg41=0.04166666666666666908*reg148; reg42=0.035220810900864524453*reg148; reg46=0.024056261216234409915*reg148; reg55=0.035220810900864524453*reg138;
   reg140=reg134*reg140; reg127=reg115*reg127; reg141=reg116-reg141; reg56=0.04166666666666666908*reg138; reg68=reg42+reg55;
   reg70=0.024056261216234409915*reg141; reg26=reg26-reg41; reg72=0.13144585576580215187*reg141; reg55=reg55+reg37; reg8=reg8-reg56;
   reg37=reg37+reg15; reg74=0.035220810900864524453*reg141; reg35=reg41+reg35; reg41=0.024056261216234395431*reg141; reg15=reg42+reg15;
   reg140=reg127-reg140; reg42=0.04166666666666666908*reg141; reg56=reg46+reg56; reg46=0.035220810900864524453*reg140; reg55=reg72+reg55;
   reg80=0.13144585576580215187*reg140; reg37=reg37+reg74; reg15=reg74+reg15; reg72=reg68+reg72; reg70=reg26-reg70;
   reg26=0.024056261216234409915*reg140; reg8=reg8-reg42; reg68=0.04166666666666666908*reg140; reg35=reg41-reg35; reg41=0.024056261216234395431*reg140;
   reg42=reg56+reg42; Ne(0,9)+=reg70-reg68; Ne(1,10)+=reg70-reg68; Ne(2,11)+=reg70-reg68; Ne(0,12)+=reg72+reg80;
   Ne(1,13)+=reg72+reg80; Ne(2,14)+=reg72+reg80; Ne(0,6)+=reg8-reg26; Ne(1,7)+=reg8-reg26; Ne(2,8)+=reg8-reg26;
   Ne(0,15)+=reg55+reg46; Ne(1,16)+=reg55+reg46; Ne(2,17)+=reg55+reg46; Ne(0,3)+=reg35-reg68; Ne(1,4)+=reg35-reg68;
   Ne(2,5)+=reg35-reg68; Ne(0,0)+=reg41-reg42; Ne(1,1)+=reg41-reg42; Ne(2,2)+=reg41-reg42; Ne(0,18)+=reg46+reg37;
   Ne(1,19)+=reg46+reg37; Ne(2,20)+=reg46+reg37; Ne(0,21)+=reg80+reg15; Ne(1,22)+=reg80+reg15; Ne(2,23)+=reg80+reg15;

}

};
} // namespace LMT

