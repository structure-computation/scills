###########################################
# Definition des variables
###########################################
#champ de deplacement par noeud
dep = Variable( unknown=True, nb_dim=[dim], default_value='0.0', unit='mm' )

#order_integration=2
integration_totale=False

#############################
# Fonctions diverses
############################

#stockage de sigma et epsilon [Sxx,Syy,Szz,Sxy,Syz,Sxz]

#assemblage du second membre et de la matrice de rigidite du maillage EF
def formulation():
  print e.interpolation["nodal"]
  res=0
  return res 

#calcul des deformations et contraintes apres avoir calcule le deplacement pour chaque noeud du maillage
def apply_on_elements_after_solve(unk_subs): # return a string
  cw = Write_code('T')
  return cw.to_string()

