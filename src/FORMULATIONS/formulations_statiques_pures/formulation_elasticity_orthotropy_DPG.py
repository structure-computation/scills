###########################################
# Definition des variables
###########################################
#forces volumiques : non utilisees pour l'instant
f_vol = Variable( interpolation='global', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )
density = Variable( interpolation='global', default_value='1', unit='kg/m^3' )
#modules d'young, coefficients de poisson et modules de cisaillement
elastic_modulus_1 = Variable( interpolation='global', default_value='157e6', unit='N/mm^2')
elastic_modulus_2 = Variable( interpolation='global', default_value='80e6', unit='N/mm^2' )
elastic_modulus_3 = Variable( interpolation='global', default_value='80e6' , unit='N/mm^2')
poisson_ratio_12 = Variable( interpolation='global', default_value='0.33' , unit='')
poisson_ratio_13 = Variable( interpolation='global', default_value='0.33' , unit='')
poisson_ratio_23 = Variable( interpolation='global', default_value='0.33' , unit='')
shear_modulus_12 = Variable( interpolation='global', default_value='50e6' , unit='N/mm^2')
shear_modulus_13 = Variable( interpolation='global', default_value='50e6' , unit='N/mm^2')
shear_modulus_23 = Variable( interpolation='global', default_value='50e6', unit='N/mm^2' )
#champ de temperature (normalement recopie du champ de temperature global sur la structure
deltaT = Variable( interpolation='global', default_value='0', unit='degC' )
epshorsplan = Variable( interpolation='global',nb_dim=[3], default_value='0.0, 0.0, 0.0', unit='' )
#coefficient de dilatation
alpha_1 = Variable( interpolation='global', default_value='0', unit='' )
alpha_2 = Variable( interpolation='global', default_value='0', unit='' )
alpha_3 = Variable( interpolation='global', default_value='0', unit='' )
#direction d'orthotropie : orthonormee prealablement
v1 = Variable( interpolation='global', nb_dim=[3], default_value='1.0, 0.0, 0.0', unit='' )
v2 = Variable( interpolation='global', nb_dim=[3], default_value='0.0, 1.0, 0.0', unit='' )

#champ de deplacement par noeud
dep = Variable( unknown=False, nb_dim=[dim], default_value='0.0', unit='mm' )
depDPG = Variable( unknown=True, nb_dim=[3], default_value='0.0', unit='mm' )
#champ de contraintes, deformations, (densite d' energie)/2 par element
sigma = Variable( interpolation='der_nodal', default_value='0', nb_dim=[6], unit='N/mm^2' )
sigma_local = Variable( interpolation='der_nodal', default_value='0', nb_dim=[6], unit='N/mm^2' )
epsilon = Variable( interpolation='der_nodal', default_value='0', nb_dim=[6], unit='' )
ener = Variable( interpolation='elementary', default_value='0', unit='N*mm' )
sigma_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[6], unit='N/mm^2' )
epsilon_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[6], unit='' )
sigma_local_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[6], unit='N/mm^2' )

#############################
# Fonctions diverses
############################

#assemblage du second membre et de la matrice de rigidite du maillage EF
def formulation():
  epsilon = grad_sym_col(depDPG.expr)
  epshorsplan3D=epshorsplan.expr
  epsilon[2]=epshorsplan3D[2]
  epsilon[4]+=epshorsplan3D[1]
  epsilon[5]+=epshorsplan3D[0]
 
  epstest = grad_sym_col(depDPG.test)

  #en 3D
  K,sigth,P =  hooke_orthotrope_th(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr,3)

  sigma = mul(K,epsilon) - sigth*deltaT.expr

  res=0
#  for i in range(1,3): res += sigma[i] * epstest3D[i]
#  i allant de 1 a 2. Par defaut, premier chiffre = 0
  for i in range(3): res += sigma[i] * epstest[i]
  for i in range(3,epsilon.size()): res += 2 * sigma[i] * epstest[i]

# energie de deformation volumique
  return res * dV

#calcul des deformations et contraintes apres avoir calcule le deplacement pour chaque noeud du maillage
def apply_on_elements_after_solve(unk_subs): # return a string
  epsilon = grad_sym_col(depDPG.expr)
  epshorsplan3D=epshorsplan.expr
  epsilon[2]=epshorsplan3D[2]
  epsilon[4]+=epshorsplan3D[1]
  epsilon[5]+=epshorsplan3D[0]
  #en 3D 
  K,sigth,P =  hooke_orthotrope_th(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr,3)

  sigma = mul(K,epsilon) - sigth*deltaT.expr

  ener=0
  for i in range(3): ener += sigma[i] * epsilon[i]
  for i in range(3,epsilon.size()): ener += 2. * sigma[i] * epsilon[i]
  ener = e.integration( ener, 2 )/2
  sigmalocal = mul(P,sigma) 
  
  #TODO dans le cas ou l'on a plusieurs points de gauss.......
  my_subs = unk_subs
  #my_subs[ time ] = time_steps[0]
  for vi in e.var_inter: my_subs[vi] = number(0.5)
  sigma = sigma.subs(EM(my_subs))
  epsilon = epsilon.subs(EM(my_subs))
  sigmalocal = sigmalocal.subs(EM(my_subs))
  
  sigma2d = vector([sigma[0],sigma[1],sigma[3]])
  eps2d = vector([epsilon[0],epsilon[1],epsilon[3]])
  sigmalocal2d = vector([sigmalocal[0],sigmalocal[1],sigmalocal[3]])
  
  if dim==2:
     cw = Write_code('T')
     for i in range(dim):
       cw.add( eps2d[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
       cw.add( sigma2d[i], 'elem.sigma[0]['+str(i)+']', Write_code.Set )
       cw.add( sigmalocal2d[i], 'elem.sigma_local[0]['+str(i)+']', Write_code.Set )
     cw.add( ener, 'elem.ener', Write_code.Set )  
     return cw.to_string()
  elif dim==3:
     cw = Write_code('T')
     for i in range(dim):
       cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
       cw.add( sigma[i], 'elem.sigma[0]['+str(i)+']', Write_code.Set )
       cw.add( sigmalocal[i], 'elem.sigma_local[0]['+str(i)+']', Write_code.Set )
     cw.add( ener, 'elem.ener', Write_code.Set )  
     return cw.to_string()
      
