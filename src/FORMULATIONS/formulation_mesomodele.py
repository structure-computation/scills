# -*- coding: utf-8 -*-
# Formulation acceptant :
# - un comportement elastique orthotrope
# - un endommagement isotrope avec comportement type diffus mesomodele
# - une deformation plastique avec ecrouissage isotrope
# - un couplage endo/plasticite via grandeurs effectives (type diffus mesomodele)


###########################################
# Definition des variables
###########################################

# ------------------------- #
# DEFINITION DES PARAMETRES #
# ------------------------- #

# ELASTICITE
# ---------
# Modules de Young
elastic_modulus_1 = Variable( interpolation='global', default_value='157e3', unit='N/mm^2')
elastic_modulus_2 = Variable( interpolation='global', default_value='8500', unit='N/mm^2' )
elastic_modulus_3 = Variable( interpolation='global', default_value='8500' , unit='N/mm^2')
# Coefficients de poisson
poisson_ratio_12 = Variable( interpolation='global', default_value='0.29' , unit='')
poisson_ratio_13 = Variable( interpolation='global', default_value='0.29' , unit='')
poisson_ratio_23 = Variable( interpolation='global', default_value='0.40' , unit='')
# Modules de cisaillement
shear_modulus_12 = Variable( interpolation='global', default_value='5000' , unit='N/mm^2')
shear_modulus_13 = Variable( interpolation='global', default_value='5000' , unit='N/mm^2')
shear_modulus_23 = Variable( interpolation='global', default_value='3000', unit='N/mm^2' )
# Masse volumique
density = Variable( interpolation='global', default_value='1', unit='kg/m^3' )

# PLASTICITE
# ----------
k_p =  Variable( interpolation='global', default_value='100e3', unit='N/mm^2')
m_p =  Variable( interpolation='global', default_value='1.0', unit='')
R0 =  Variable( interpolation='global', default_value='0', unit='N/mm^2')
couplage = Variable( interpolation='global', default_value='0', unit='')

# ENDOMMAGEMENT
# -------------
Yo = Variable( interpolation='global', default_value="0.0", unit='N/mm^2')
Yc = Variable( interpolation='global', default_value="1.0", unit='N/mm^2')
Ycf = Variable( interpolation='global', default_value="1.0", unit='N/mm^2')
dmax = Variable( interpolation='global', default_value="0.999", unit='')
a = Variable( interpolation='global', default_value="1.0", unit='')
tau_c = Variable( interpolation='global', default_value="1.0", unit='')
b_c = Variable( interpolation='global', default_value="1.0", unit='')
effet_retard = Variable( interpolation='global', default_value="0", unit='')
pli = Variable ( interpolation='elementary', default_value='0', unit='')

# THERMIQUE
# ---------
# Coefficients de dilatation
alpha_1 = Variable( interpolation='global', default_value='2.8e-6', unit='' )
alpha_2 = Variable( interpolation='global', default_value='30e-6', unit='' )
alpha_3 = Variable( interpolation='global', default_value='30e-6', unit='' )

# CALCUL
# ------
# Directions d'orthotropie
v1 = Variable( interpolation='global', nb_dim=[3], default_value='1.0, 0.0, 0.0', unit='' )
v2 = Variable( interpolation='global', nb_dim=[3], default_value='0.0, 1.0, 0.0', unit='' )
# Contrainte (1) ou Deformation (0) plane
resolution = Variable( interpolation='global', default_value='0', unit='' )
# Integration totale (?)
integration_totale=False


# ------------------------ #
# DEFINITION DES VARIABLES #
# ------------------------ #

#forces volumiques : non utilisees pour l'instant
f_vol = Variable( interpolation='global', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )

#champ de temperature (normalement recopie du champ de temperature global sur la structure
deltaT = Variable( interpolation='global', default_value='0', unit='degC' )


#champ de deplacement par noeud
dep = Variable( unknown=True, nb_dim=[dim], default_value='0.0', unit='mm' )

#champ de contraintes, deformations, (densite d' energie)/2 par element
sigma = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_local = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
epsilon = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
ener = Variable( interpolation='elementary', default_value='0', unit='N*mm' )
sigma_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
epsilon_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
sigma_local_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_von_mises = Variable( interpolation='elementary', default_value='0',unit='N/mm^2' )

# Variables de plasticite
epsilon_e = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
epsilon_p = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
p = Variable( interpolation='elementary', default_value='0', unit='' )
R_p = Variable( interpolation='elementary', default_value='0', unit='N/mm^2' )

# Variables d'endommagement
d = Variable( interpolation='elementary', default_value='0', nb_dim=[3], unit='' )
Yd = Variable( interpolation='elementary', default_value='0', nb_dim=[3], unit='N/mm^2' )
calcul_rigi = k_p =  Variable( interpolation='global', default_value='0', unit='')

#############################
# Fonctions diverses
############################

#stockage de sigma et epsilon [Sxx,Syy,Szz,Sxy,Syz,Sxz]

#assemblage du second membre et de la matrice de rigidite du maillage EF
def formulation():
  epsilon = grad_sym_col(dep.expr)
  epstest = grad_sym_col(dep.test)
  #en 3D
  endo = d.expr
  endo0 = endo[0]
  endo1 = endo[1]
  endo2 = endo[2]
  K0,H0,epsth0,P0 =  hooke_orthotrope_th_3d_endo(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr,endo0,endo1,endo2)
  
  if (dim==2):
      #simplifications dimension ou comportement
      K0s,epsth0s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane stress')
      K1s,epsth1s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane strain')
      #P=simplification_projection(P0,dim)
      
      K= K0s*resolution.expr + K1s*(1-resolution.expr)
      epsth= epsth0s*resolution.expr + epsth1s*(1-resolution.expr)
  elif(dim==3):
      K= K0
      epsth = epsth0
  sigma = mul(K, ( epsilon - epsilon_p.expr - epsth*deltaT.expr) )
  
  #sm_p = mul(K, ( epsilon_p.expr ) )
  
  res=0
  for i in range(dim): res += sigma[i] * epstest[i]
  for i in range(dim,epsilon.size()): res += 2 * sigma[i] * epstest[i]
  
  #sm = 0
  #for i in range(dim): sm += sm_p[i] * epstest[i]
  #for i in range(dim,epsilon.size()): sm += 2 * sm_p[i] * epstest[i]

  return (res) * dV


#calcul des deformations et contraintes apres avoir calcule le deplacement pour chaque noeud du maillage
def apply_on_elements_after_solve(unk_subs): # return a string
  epsilon = grad_sym_col(dep.expr)
  epsilon_e = epsilon - epsilon_p.expr
  #en 3D
  endo = d.expr
  endo0 = endo[0]
  endo1 = endo[1]
  endo2 = endo[2]
  K0,H0,epsth0,P0 =  hooke_orthotrope_th_3d_endo(elastic_modulus_1.expr, elastic_modulus_2.expr, elastic_modulus_3.expr, poisson_ratio_12.expr, poisson_ratio_13.expr, poisson_ratio_23.expr, shear_modulus_12.expr, shear_modulus_13.expr, shear_modulus_23.expr,  v1.expr, v2.expr, alpha_1.expr, alpha_2.expr, alpha_3.expr, 0, 0, 0)
  
  if (dim==2):
      #simplifications dimension ou comportement
      K0s,epsth0s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane stress')
      K1s,epsth1s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane strain')
      P=simplification_projection(P0,dim)
      
      K= K0s*resolution.expr + K1s*(1-resolution.expr)
      epsth= epsth0s*resolution.expr + epsth1s*(1-resolution.expr)
  elif (dim==3):
      K = K0
      epsth = epsth0
      P = P0
  
  sigma = mul(K, ( epsilon_e - epsth*deltaT.expr) )
  
  ener=0
  for i in range(dim): ener += sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
  for i in range(dim,epsilon.size()): ener += 2 * sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
  ener = e.integration( ener, 2 )/2.
  
  sigmalocal = mul(P,sigma) 
  sigma_3d,epsilon_3d = reconstruction_quantites_3d(epsilon,K0,H0,epsth0,deltaT.expr,dim,type_stress_2D='plane stress')
  sigma_von_mises=von_mises( sigma_3d )
  
  #TODO dans le cas ou l'on a plusieurs points de gauss.......
  my_subs = unk_subs
  #my_subs[ time ] = time_steps[0]
  for vi in e.var_inter: my_subs[vi] = number(0.5)
  sigma = sigma.subs(EM(my_subs))
  epsilon = epsilon.subs(EM(my_subs))
  epsilon_e = epsilon_e.subs(EM(my_subs))
  sigmalocal = sigmalocal.subs(EM(my_subs))
  sigma_von_mises = sigma_von_mises.subs(EM(my_subs))
  
  cw = Write_code('T')
  for i in range(dim*(dim+1)/2):
    cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
    cw.add( epsilon_e[i], 'elem.epsilon_e[0]['+str(i)+']', Write_code.Set )
    #cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
    cw.add( sigma[i], 'elem.sigma[0]['+str(i)+']', Write_code.Set )
    cw.add( sigmalocal[i], 'elem.sigma_local[0]['+str(i)+']', Write_code.Set )
  cw.add( ener, 'elem.ener', Write_code.Set )  
  cw.add( sigma_von_mises, 'elem.sigma_von_mises', Write_code.Set )  
  return cw.to_string()

#calcul des deformations et contraintes uniquement
def apply_on_elements_after_solve_2(unk_subs): # return a string
  epsilon = grad_sym_col(dep.expr)
  epsilon_e = epsilon - epsilon_p.expr
  #en 3D
  endo = d.expr
  endo0 = endo[0]
  endo1 = endo[1]
  endo2 = endo[2]
  K0,H0,epsth0,P0 =  hooke_orthotrope_th_3d_endo(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr,0,0,0)
  
  if (dim==2):
      #simplifications dimension ou comportement
      K0s,epsth0s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane stress')
      K1s,epsth1s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane strain')
      P=simplification_projection(P0,dim)
      
      K= K0s*resolution.expr + K1s*(1-resolution.expr)
      epsth= epsth0s*resolution.expr + epsth1s*(1-resolution.expr)
  elif (dim==3):
      K = K0
      epsth = epsth0
      P = P0
  
  sigma = mul(K, ( epsilon_e - epsth*deltaT.expr) )
  
  sigmalocal = mul(P,sigma) 
  sigma_3d,epsilon_3d = reconstruction_quantites_3d(epsilon,K0,H0,epsth0,deltaT.expr,dim,type_stress_2D='plane stress')
  sigma_von_mises=von_mises( sigma_3d )
  
  #TODO dans le cas ou l'on a plusieurs points de gauss.......
  my_subs = unk_subs
  #my_subs[ time ] = time_steps[0]
  for vi in e.var_inter: my_subs[vi] = number(0.5)
  sigma = sigma.subs(EM(my_subs))
  epsilon = epsilon.subs(EM(my_subs))
  epsilon_e = epsilon_e.subs(EM(my_subs))
  
  cw = Write_code('T')
  for i in range(dim*(dim+1)/2):
    cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
    cw.add( epsilon_e[i], 'elem.epsilon_e[0]['+str(i)+']', Write_code.Set )
    #cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
    cw.add( sigma[i], 'elem.sigma[0]['+str(i)+']', Write_code.Set )
  return cw.to_string()