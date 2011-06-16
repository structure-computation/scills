Import('env pb_libs ')

liste1 = [ 
   './src/multiscale.cpp',
   './src/MAILLAGE/multiscale_geometry_mesh.cpp',
   './src/MATERIAUX/assignation_materials_property.cpp',
   './src/OPERATEURS/multiscale_operateurs.cpp',
   './src/ITERATIONS/iterate_stat_Qstat.cpp',
   './src/POSTTRAITEMENTS/affichage.cpp',
]

liste2 = [
#    './ITERATIONS/iterate_stat_incrementale.cpp'
#       './PROBMICRO/calcul_cohesif.cpp',
#    './PROBMICRO/preselection_interface.cpp',
#    './PROBMICRO/selection_interface.cpp',
]

listetot = liste1 + liste2

env.Depends( listetot, pb_libs )

liste = env.Library('mon_prog',listetot)

Return('liste')
