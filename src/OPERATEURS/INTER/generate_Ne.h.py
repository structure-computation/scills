from LMT.formal_lf import *
import string
import sys

nbelem={}
nbelem[0]=Element('Bar',2)
nbelem[1]=Element('Triangle',3)
nbelem[2]=Element('Quad',3)
nbelem[3]=Element('Bar_3',2)
nbelem[4]=Element('Triangle_6',3)
nbelem[5]=Element('Quad_8',3)

print """
#include "mesh/bar.h"
#include "mesh/triangle.h"
#include "mesh/quad.h"
#include "mesh/bar_3.h"
#include "mesh/triangle_6.h"
#include "mesh/quad_8.h"

namespace LMT {
struct add_elem_Ne{
"""

for l,integ in zip(range(0,6),[0,1,2,1,2,3]):
   e =  nbelem[l]
   
   # calcul des termes de la matrices de sous integration
   # Nij=int(phi_i(X)dV)
   
   fctbase={}
   for i in range(0,e.nb_nodes):
      fctbase[i]=e.interpolation['nodal'].diff(e.val[i])
   
   N=ExVector(e.nb_nodes)
      
   for i in range(0,e.nb_nodes):
      N[i] = e.integration(fctbase[i],integ) 
   
   # ecriture du code
   cw = Write_code('double',3)
   
   print """ /// pour les %s
   template<class TNB,class TN,class TD,unsigned NET,class TMAT>
   void operator()(const Element<%s,TNB,TN,TD,NET> e,TMAT &Ne) const {
   """%(e.name,e.name)
   
   for j in range(0,e.dim):
      for i in range(0,e.nb_nodes):
         tot=N[i]
         cw.add( tot, 'Ne(%i,%i)'%(j,j+i*e.dim), Write_code.Add )
   
   print string.replace( cw.to_string(), 'elem.', 'e.' )
   
   print "}"
print """
};
} // namespace LMT
"""
