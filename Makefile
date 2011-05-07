all:local 

codegen_py:
	cd LMT/include/codegen; scons

local: codegen_py
	scons -j4 dep_py=1 arch=x86-64

clean:
	scons -c
	cd LMT/include/codegen; scons -c

vouvray: codegen_py
	rsync --exclude '*~' --exclude '*.o' --exclude 'tmp' --exclude 'build' --exclude 'multinew*'  -auv -e ssh . caignot@vouvray:/home/caignot/MULTI_MPI
	ssh caignot@vouvray make

vouvrayviamadiran:
	rsync --exclude '*~' --exclude '*.o' --exclude 'tmp' --exclude 'LMT.old' --exclude '_darcs/' --exclude 'build' --exclude 'multinew*' --include '*.py'  -auv -e ssh . caignot@ssh.lmt.ens-cachan.fr:/u/caignot/_CPP/MULTI_v09
	ssh caignot@ssh.lmt.ens-cachan.fr make

p4viamadiran:
	cd src;rsync --exclude '*~' --exclude '*.o' --exclude 'tmp' --exclude 'multinew*'  -auv -e ssh . caignot@ssh.lmt.ens-cachan.fr:/u/caignot/_CPP/MULTI_MPI/src
	ssh caignot@ssh.lmt.ens-cachan.fr 'make p4'

kurofu: codegen_py
	rsync --exclude '*~' --exclude '*.o' --exclude 'mat' -au -e ssh . leclerc@kurofu.meso.ens-cachan.fr:/home/people/leclerc/fm
	ssh leclerc@kurofu.meso.ens-cachan.fr ". .bash_profile; cd fm; scons -j4 dep_py=0 debug=0 opt=1 arch=ia64"

maekake: codegen_py
	rsync --exclude '*~' --exclude '*.o' -au -e ssh . caignot@maekake.meso.ens-cachan.fr:/home/people/caignot/MULTI

local1:
	g++ treillis.cpp -ILMT/include -O3 -o treillis

rsync:
	rsync -auv -e ssh . alain@alain201180.no-ip.info:/home/alain/LMT/_CPP/MULTI_v09_MPI

rsync_lmt:
	rsync --exclude 'tmp' -auv -e ssh ./src alain@portalain:/home/alain/LMT/_CPP/MULTI_v09_MPI
	rsync --exclude 'tmp' -auv -e ssh ./build alain@portalain:/home/alain/LMT/_CPP/MULTI_v09_MPI

rsync_madiran:
	rsync --exclude 'tmp' -auv -e ssh . caignot@ssh.lmt.ens-cachan.fr:/u/caignot/_CPP/MULTI_v09

rsync_home:
	rsync --exclude 'tmp' --exclude '*~' -auv -e ssh ./src alain@alain201180.no-ip.info:/home/alain/LMT/_CPP/MULTI_v09_MPI
	rsync --exclude 'tmp' --exclude '*~' -auv -e ssh ./build alain@alain201180.no-ip.info:/home/alain/LMT/_CPP/MULTI_v09_MPI

