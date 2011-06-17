##############################################
# Procedure pour calculer la mesure d'un element
##############################################

from LMT.formal_lf import *
import string
import sys

# choix de l'element
e=Element('Hexa',3)
print """
namespace LMT {
template<class TN,class TNG,class TD,unsigned NET>
typename TNG::T measure( const Element<%s,TN,TNG,TD,NET> &elem ) {
typedef typename TNG::T T;
"""%e.name
print " T V=0.;"

#calcul de la mesure
V = e.integration(number(1.),2)

cw = Write_code('T',3)
cw.add( V, 'V', Write_code.Add )
print string.replace( cw.to_string(), 'V', 'V' )
print """
return V;
}
};
"""
    