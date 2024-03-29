###########################################
# Definition des variables
###########################################
#forces volumiques : non utilisees pour l'instant
f_vol = Variable( interpolation='global', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )
f_vol_e = Variable( interpolation='elementary', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )
density = Variable( interpolation='global', default_value='1', unit='kg/m^3' )
#modules d'young, coefficients de poisson et modules de cisaillement
elastic_modulus_1 = Variable( interpolation='global', default_value='157e3', unit='N/mm^2')
elastic_modulus_2 = Variable( interpolation='global', default_value='8500', unit='N/mm^2' )
elastic_modulus_3 = Variable( interpolation='global', default_value='8500' , unit='N/mm^2')
poisson_ratio_12 = Variable( interpolation='global', default_value='0.29' , unit='')
poisson_ratio_13 = Variable( interpolation='global', default_value='0.29' , unit='')
poisson_ratio_23 = Variable( interpolation='global', default_value='0.40' , unit='')
shear_modulus_12 = Variable( interpolation='global', default_value='5000' , unit='N/mm^2')
shear_modulus_13 = Variable( interpolation='global', default_value='5000' , unit='N/mm^2')
shear_modulus_23 = Variable( interpolation='global', default_value='3000', unit='N/mm^2' )
#champ de temperature (normalement recopie du champ de temperature global sur la structure
#deltaT = Variable( interpolation='global', default_value='0', unit='degC' )
deltaT = Variable( interpolation='elementary', default_value='0.0', unit='degC' )
#coefficient de dilatation
alpha_1 = Variable( interpolation='global', default_value='2.8e-6', unit='' )
alpha_2 = Variable( interpolation='global', default_value='30e-6', unit='' )
alpha_3 = Variable( interpolation='global', default_value='30e-6', unit='' )
#direction d'orthotropie : orthonormee prealablement
#angle = Variable( interpolation='elementary', default_value='0', unit='' )
v1 = Variable( interpolation='global', nb_dim=[3], default_value='1.0, 0.0, 0.0', unit='' )
v2 = Variable( interpolation='global', nb_dim=[3], default_value='0.0, 1.0, 0.0', unit='' )

#booleen valant 1 pour une formulation en contrainte plane 0 pour deformation plane
resolution = Variable( interpolation='global', default_value='0', unit='' )

#champ de deplacement par noeud
dep = Variable( unknown=True, nb_dim=[dim], default_value='0.0', unit='mm' )

#champ de contraintes, deformations, (densite d' energie)/2 par element
sigma = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_local = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
epsilon = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
ener = Variable( interpolation='elementary', default_value='0', unit='N*mm' )
sigma_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
epsilon_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
sigma_mises_skin = Variable( interpolation='skin_elementary', default_value='0',  unit='N/mm^2' )
sigma_local_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_von_mises = Variable( interpolation='elementary', default_value='0',unit='N/mm^2' )
sigma_mises_skin = Variable( interpolation='skin_elementary', default_value='0', unit='N/mm^2' )

integration_totale=False

#############################
# Fonctions diverses
############################

#stockage de sigma et epsilon [Sxx,Syy,Szz,Sxy,Syz,Sxz]

#assemblage du second membre et de la matrice de rigidite du maillage EF
def formulation():
  epsilon = grad_sym_col(dep.expr)
  epstest = grad_sym_col(dep.test)
  #en 3D
  K0,H0,epsth0,P0 =  hooke_orthotrope_th_3d(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr)
  
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
  
  sigma = mul(K, ( epsilon - epsth*deltaT.expr) )
  
  res=0
  for i in range(dim): res += sigma[i] * epstest[i]
  for i in range(dim,epsilon.size()): res += 2 * sigma[i] * epstest[i]
  res-=dot( f_vol_e.expr , dep.test )
  
  #res -= dot( f_vol.expr , dep.test )

  return res * dV


#calcul des deformations et contraintes apres avoir calcule le deplacement pour chaque noeud du maillage
def apply_on_elements_after_solve(unk_subs): # return a string
  epsilon = grad_sym_col(dep.expr)
  #en 3D
  K0,H0,epsth0,P0 =  hooke_orthotrope_th_3d(elastic_modulus_1.expr,elastic_modulus_2.expr,elastic_modulus_3.expr,
  poisson_ratio_12.expr,poisson_ratio_13.expr,poisson_ratio_23.expr,
  shear_modulus_12.expr,shear_modulus_13.expr,shear_modulus_23.expr,
  v1.expr,v2.expr,alpha_1.expr,alpha_2.expr,alpha_3.expr)
  
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
  
  sigma = mul(K, ( epsilon - epsth*deltaT.expr) )
  
  ener=0
  for i in range(dim): ener += sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
  for i in range(dim,epsilon.size()): ener += 2 * sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
  ener = e.integration( ener, 2 )/2.
  
  sigmalocal = mul(P,sigma) 
  epsilon_3d, sigma_3d = reconstruction_quantites_3d(epsilon,K0,H0,epsth0,deltaT.expr,dim,type_stress_2D='plane stress')
  sigma_von_mises=von_mises( sigma_3d )
  
  #TODO dans le cas ou l'on a plusieurs points de gauss.......
  my_subs = unk_subs
  #my_subs[ time ] = time_steps[0]
  for vi in e.var_inter: my_subs[vi] = number(0.5)
  sigma = sigma.subs(EM(my_subs))
  epsilon = epsilon.subs(EM(my_subs))
  sigmalocal = sigmalocal.subs(EM(my_subs))
  sigma_von_mises = sigma_von_mises.subs(EM(my_subs))
  
  cw = Write_code('T')
  for i in range(dim*(dim+1)/2):
    cw.add( epsilon[i], 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
    cw.add( sigma[i], 'elem.sigma[0]['+str(i)+']', Write_code.Set )
    cw.add( sigmalocal[i], 'elem.sigma_local[0]['+str(i)+']', Write_code.Set )
  cw.add( ener, 'elem.ener', Write_code.Set )  
  cw.add( sigma_von_mises, 'elem.sigma_von_mises', Write_code.Set )  
  return cw.to_string()

