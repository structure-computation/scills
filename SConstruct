#from fetch import *
#if not os.access('LMT',os.F_OK): fetch_LMT('','','')
from LMT import *
from string import *

#opt = 'opt' in ARGUMENTS and int(ARGUMENTS['opt'])
#debug = 'debug' in ARGUMENTS and int(ARGUMENTS['debug'])
#timdavis = 'timdavis' in ARGUMENTS and int(ARGUMENTS['timdavis'])
if 'arch' in ARGUMENTS: arch = ARGUMENTS['arch']
else: arch=''

#########################################################################
# Parametres a changer
#########################################################################

#Mettre le flag opt a 1 pour compiler avec optimisation (plus lent) sinon en mettant a 0 on compile en mode debug
opt = 1 
debug = 1-opt
timdavis=0 #pour utiliser un solveur different : pas encore bien teste

#Repertoires contenant les .cpp et .h
mes_repertoires = 'DEFINITIONS  POSTTRAITEMENTS  OPERATEURS/SST  OPERATEURS/INTER  OPERATEURS/MACRO MAILLAGE  MATERIAUX ITERATIONS ITERATIONS/LINEAR  ITERATIONS/LOCAL  ITERATIONS/ERROR   PROBMICRO  UTILITAIRES  FORMULATIONS  ALLOCATIONS    MPI GEOMETRY  COMPUTE '
mes_repertoires = map( lambda x: 'src/'+x, split(mes_repertoires))

list_repertoires = ['LMT/include', 'build']+mes_repertoires + ['LMT/formulations']
list_repertoires_cppflags = map( lambda x: os.getcwd()+'/'+x, list_repertoires)

#choix des formulations a compiler
formuortho = 1

formuiso   = 1
formuisop  = 1
formuendo  = 1
formumeso  = 1


#pour avoir les cout dans un ficher ...
crout=1

#pour compiler uniquement en 2 ou 3d (mettre les 2 a 1 pour compiler en 2D et en 3D
dim2 = 0
dim3 = 1

mes_formulations =  ['elasticity_isotropy_stat_Qstat']*formuiso +  ['elasticity_orthotropy_stat_Qstat']*formuortho + ['plasticity_isotropy_stat_Qstat']*formuisop + ['mesomodele']*formumeso + ['elasticity_damageable_isotropy_stat_Qstat']*formuendo

#pour utiliser la p-surdiscretisation mettre le flag a 1 (on compile alors avec les elements indiques + les elements de degre 2)         
sur_discretisation=0
#choix des elements ('Triangle', 'Quad', 'Tetra', 'Wedge', 'Hexa')
mes_elements = [ 'Triangle', 'Quad', 'Tetra', 'Wedge', 'Hexa']

# Les fichiers .o sont places dans build/ et dans le repertoire debug/opt ainsi que dans monprog_dim2ou3
nom_repertoire_prog='build/'+arch+'_debug'*debug+'_opt'*opt+'/monprog_dim'+'2'*dim2+'_'+dim3*'3'

# generation automatique des certains fichiers du programme en fonction des formulations
# Obsolete
#import src.generation_auto
#src.generation_auto.generation_auto(mes_formulations)
#######################################################################



if(sur_discretisation==1):
   for i in range(len(mes_elements)):
      if(mes_elements[i]=="Quad"):
         mes_elements+=["Quad_8"]
      elif(mes_elements[i]=="Triangle"):
         mes_elements+=["Triangle_6"]
      elif(mes_elements[i]=="Tetra"):
         mes_elements+=["Tetra_10"]
      elif(mes_elements[i]=="Wedge"):
         mes_elements+=["Wedge_15"]
      elif(mes_elements[i]=="Hexa"):
         mes_elements+=["Hexa_20"]

print mes_elements


# --------------------------------------------------------------------------------------------


# compiler choice
if os.access('/opt/gcc-4.0/bin/g++',os.F_OK): comp,compp = '/opt/gcc-4.0/bin/g++', '/opt/gcc-4.0'
else: comp,compp = 'g++', '/usr'

if 'icc' in comp or 'icpc' in comp:
    cppf = ' -g ' + \
           ' -DNOx86 ' * (arch=='ia64') + \
           ' -O3 ' * opt + \
           ' -g -DDEBUG ' * debug + \
           (' -mcpu='+arch+' ')*( len(arch)!=0 ) # 
else:
    cppf = ' -pipe -Wall -Wno-deprecated -w' + \
           ' -DNOx86 ' * (arch=='ia64') + \
           ' -O3 ' * opt + \
           ' -ffast-math -funroll-loops -fexpensive-optimizations ' * ( opt and arch!='ia64' ) + \
           ' -g3 -DDEBUG ' * debug + \
           (' -march='+arch+' ')*( len(arch)!=0 and arch!='ia64' ) # 


if (timdavis==1):
  cppftd = file("/opt/MATH/timdavis/cppflags","r").read()
  ldflagtd = file("/opt/MATH/timdavis/ldflags","r").read()
else:
  cppftd = ''
  ldflagtd = ''

env = Environment(
    CXX = comp,
    CPPPATH = list_repertoires_cppflags,
#    LINKFLAGS = ' -LUTIL/lam/lib -llammpi++ -lmpi -llam -lutil -ldl -L'+compp+'/lib'+'64'*(arch=='x86-64')+'/ ' + linkflags( ['xml2-config'] )+ldflagtd,
    LINKFLAGS = ' -LUTIL/lam/lib -lmpi_cxx -lmpi -lutil -ldl -L'+compp+'/lib'+'64'*(arch=='x86-64')+'/ ' + linkflags( ['xml2-config'] )+ldflagtd,
    CPPFLAGS = cppf + cppflags( ['xml2-config'] )+' -DLDL '*(1-timdavis)+cppftd+' -DDIMENSION3'*dim3+' -DDIMENSION2'*dim2+' -DFORMUORTHO'*formuortho+'  -Dcrout_alain'*crout+' -IUTIL/openmpi ',
)
# -L/usr/lib/lam/lib -llammpi++ -lmpi -llam 
# -LUTIL/lam/lib -llammpi++ -lmpi -llam 
env.Command('fetch','',fetch_LMT)
env.Command('dist','',makedist)
if ( not 'dep_py' in ARGUMENTS ) or int(ARGUMENTS['dep_py']):
    make_dep_py(env)

# from LMT.formal_lf.variable import *

from LMT import *
from LMT.formal_lf import *

# Modifiable
# formulations et elements modifiables, possibilite d'ajouter des champs par Sst
# on doit pouvoir passer les vrai type en faisant T='int' ...
pb_libs=make_pb(
    env,
    opt,
    name = 'pb_elast',
    incpaths=list_repertoires_cppflags,
    formulations = mes_formulations ,
    elements = mes_elements,
    additional_fields = {
        "qtrans" : Variable(unit='mm',T='Vec<double,DIM>', default_value='0'),
        "typmat" : Variable(interpolation='elementary',unit='',default_value='0'),
        "numsst" : Variable(interpolation='elementary',unit='',default_value='0'),
        "num_proc" : Variable(interpolation='elementary',unit='',default_value='0'),
        "typmat_skin" : Variable(interpolation='skin_elementary',unit='',default_value='0'),
        "numsst_skin" : Variable(interpolation='skin_elementary',unit='',default_value='0'),
        "num_proc_skin" : Variable(interpolation='skin_elementary',unit='',default_value='0'),
},
   dep_py = 1
)

BuildDir(nom_repertoire_prog+'/LMT', 'LMT/include' ,duplicate=0)
libs = SConscript('LMT/include/SConscript', exports='env', build_dir=nom_repertoire_prog+'/LMT')
BuildDir(nom_repertoire_prog+'/MULTI', './' ,duplicate=0)
liste = SConscript('SConscript', exports='env pb_libs', build_dir=nom_repertoire_prog+'/MULTI')

libs += [
    '/opt/MATH/CHOLAMD/CHOLMOD/Lib/libcholmod.a',
    '/opt/MATH/CHOLAMD/UMFPACK/Lib/libumfpack.a',
    '/opt/MATH/CHOLAMD/AMD/Lib/libamd.a',
    '/opt/MATH/CHOLAMD/COLAMD/libcolamd.a',
    '/opt/MATH/CHOLAMD/CCOLAMD/libccolamd.a',
    '/opt/MATH/metis-4.0/libmetis.a',
    '/opt/MATH/CHOLAMD/UFconfig/xerbla/libcerbla.a',
]*timdavis

libs += [
    'UTIL/metis/libmetis.a',
]*(1-timdavis)

#OLD VOUVRAY
#libs += [
#     'UTIL/lam/lib/liblammpi++.a',
#     'UTIL/lam/lib/libmpi.a',
#     'UTIL/lam/lib/liblam.a',
#     '/usr/lib64/libpthread.a',
#]

#libs += [
#  'lammpi++',
#  'mpi',
#  'lam',
#  'pthread',
#]

flag=0 #0 pour compiler tous les fichiers cpp , 1 pour le fichier test seulement
autres=[
 #'test.cpp',
 #'makerefinement.cpp'
 #  'FICHIERS_TEST/conversion_maillage_data_gmsh.cpp'
# 'validation/validation_elements.cpp',
#  'test_structure.cpp',
#  'treillis.cpp',
]
autres+=libs

#pour changer la dimension

env2=env
env2["CPPFLAGS"]+=' -DTYPEREEL=double'
#prg = env2.Program(('multinew')*(1-flag)+flag*'test_struct',(pb_libs+liste+libs)*(1-flag)+flag*(pb_libs+autres))
#prg = env2.Program('compil_init',(pb_libs))
#Default(prg)
