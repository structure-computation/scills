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
elastic_modulus = Variable( interpolation='global', default_value='200e3', unit='N/mm^2')
# Coefficients de poisson
poisson_ratio = Variable( interpolation='global', default_value='0.3' , unit='')
# Masse volumique
density = Variable( interpolation='global', default_value='1', unit='kg/m^3' )

# PLASTICITE
# ----------
k_p =  Variable( interpolation='global', default_value='100e3', unit='N/mm^2')
m_p =  Variable( interpolation='global', default_value='1.0', unit='')
R0 =  Variable( interpolation='global', default_value='0', unit='N/mm^2')

# THERMIQUE
# ---------
# Coefficients de dilatation
alpha = Variable( interpolation='global', default_value='0.0', unit='' )

# CALCUL
# ------
# Contrainte (1) ou Deformation (0) plane
resolution = Variable( interpolation='global', default_value='0', unit='' )
# Integration totale (?)
integration_totale=False


# ------------------------ #
# DEFINITION DES VARIABLES #
# ------------------------ #

#forces volumiques : non utilisees pour l'instant
f_vol = Variable( interpolation='global', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )
f_vol_e = Variable( interpolation='elementary', nb_dim=[dim], default_value='0.0,'*(dim-1)+'0.0', unit='N/m^3' )

#champ de temperature (normalement recopie du champ de temperature global sur la structure
deltaT = Variable( interpolation='global', default_value='0', unit='degC' )


#champ de deplacement par noeud
dep = Variable( unknown=True, nb_dim=[dim], default_value='0.0', unit='mm' )

# Variables d'elasticite
sigma = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_local = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_local_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='N/mm^2' )
sigma_von_mises = Variable( interpolation='elementary', default_value='0',unit='N/mm^2' )
epsilon = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
epsilon_skin = Variable( interpolation='skin_elementary', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )

# Variables de plasticite
epsilon_e = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
epsilon_p = Variable( interpolation='der_nodal', default_value='0', nb_dim=[dim*(dim+1)/2], unit='' )
p = Variable( interpolation='elementary', default_value='0', unit='' )
R_p = Variable( interpolation='elementary', default_value='0', unit='N/mm^2' )

# (Densite d'energie)/2 par element
ener = Variable( interpolation='elementary', default_value='0', unit='N*mm' )

#############################
# Fonctions diverses
############################

#stockage de sigma et epsilon [Sxx,Syy,Szz,Sxy,Syz,Sxz]

#assemblage du second membre et de la matrice de rigidite du maillage EF
def formulation():
    epsilon = grad_sym_col(dep.expr)
    epstest = grad_sym_col(dep.test)
    #en 3D
    K0 ,H0, epsth0 = hooke_isotrope_th_3d(elastic_modulus.expr,poisson_ratio.expr,alpha.expr)

    if(dim==2):
        #simplification dimension ou comportement
        K0s,epsth0s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane stress')
        K1s,epsth1s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane strain')
        K= K0s*resolution.expr + K1s*(1-resolution.expr)
        epsth = epsth0s*resolution.expr + epsth1s*(1-resolution.expr)
    elif(dim==3):
        K = K0
        epsth = epsth0
        
    sigma = mul(K, (epsilon - epsilon_p.expr - epsth*deltaT.expr) )
    res=0
    for i in range(dim): res += sigma[i] * epstest[i] 
    for i in range(dim,epstest.size()): res += 2 * sigma[i] * epstest[i]
    res -= dot( f_vol_e.expr , dep.test )
    res -= dot( f_vol.expr   , dep.test )
    return res * dV


#calcul des deformations et contraintes apres avoir calcule le deplacement pour chaque noeud du maillage
def apply_on_elements_after_solve(unk_subs): # return a string
    epsilon = grad_sym_col(dep.expr)
    epsilon_e = epsilon - epsilon_p.expr
    #en 3D
    K0 ,H0, epsth0 = hooke_isotrope_th_3d(elastic_modulus.expr,poisson_ratio.expr,alpha.expr)

    if(dim==2):
        #simplification dimension ou comportement
        K0s,epsth0s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane stress')
        K1s,epsth1s = simplification_behaviour(K0,H0,epsth0,dim,type_stress_2D='plane strain')
        K= K0s*resolution.expr + K1s*(1-resolution.expr)
        epsth = epsth0s*resolution.expr + epsth1s*(1-resolution.expr)
    elif(dim==3):
        K = K0
        epsth = epsth0
        
    sigma = mul(K, (epsilon_e - epsth*deltaT.expr) )
  
    ener=0
    for i in range(dim): ener += sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
    for i in range(dim,epsilon.size()): ener += 2 * sigma[i] * (epsilon[i] - epsth[i]*deltaT.expr)
    ener = e.integration( ener, 2 )/2.
    
    sigmalocal = sigma
    epsilon_3d,sigma_3d = reconstruction_quantites_3d(epsilon,K0,H0,epsth0,deltaT.expr,dim,type_stress_2D='plane stress')
    sigma_von_mises=von_mises( sigma_3d )
    R_p = R0.expr + k_p.expr*p.expr**m_p.expr
  
    #TODO dans le cas ou l'on a plusieurs points de gauss.......
    my_subs = unk_subs
    #my_subs[ time ] = time_steps[0]
    for vi in e.var_inter: my_subs[vi] = number(0.5)
    sigma = sigma.subs(EM(my_subs))
    epsilon = epsilon.subs(EM(my_subs))
    epsilon_e = epsilon_e.subs(EM(my_subs))
    epsilon_p = epsilon_p.subs(EM(my_subs))
    sigmalocal = sigmalocal.subs(EM(my_subs))
    sigma_von_mises = sigma_von_mises.subs(EM(my_subs))
    R_p = R_p.subs(EM(my_subs))
    
    cw = Write_code('T')
    for i in range(dim*(dim+1)/2):
        cw.add( epsilon[i]   , 'elem.epsilon[0]['+str(i)+']', Write_code.Set )
        cw.add( epsilon_e[i] , 'elem.epsilon_e[0]['+str(i)+']', Write_code.Set )
        cw.add( epsilon_p[i] , 'elem.epsilon_p[0]['+str(i)+']', Write_code.Set )
        cw.add( sigma[i]     , 'elem.sigma[0]['+str(i)+']', Write_code.Set )
        cw.add( sigmalocal[i], 'elem.sigma_local[0]['+str(i)+']', Write_code.Set )
    cw.add( ener            , 'elem.ener', Write_code.Set )  
    cw.add( sigma_von_mises , 'elem.sigma_von_mises', Write_code.Set )  
    cw.add( R_p             , 'elem.R_p', Write_code.Set )  
    return cw.to_string()
