# petsc and slepc settings
SLEPC_DIR = ${HOME}/local/src/slepc-3.4.3
PETSC_DIR = ${HOME}/local/src/petsc-3.4.3
PETSC_ARCH = rccls_C

CXXFLAGS=-O3 -Wall
#$CXXFLAGS=-g

pep8:
	pep8 bspline.py

bspline_bind_linux: bspline_bind.cpp bspline.cpp bspline.hpp
	g++ -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -o bspline_bind.so bspline_bind.cpp bspline.cpp -lboost_python -lboost_numpy -L${LIBRARY_PATH}

bspline_bind_linux2: bspline_bind.cpp bspline.cpp bspline.hpp
	g++ -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -o bspline_bind.so bspline_bind.cpp bspline.cpp -lboost_python -lboost_numpy ${CXXFLAGS}


bspline_bind_mac.so: bspline_bind.cpp bspline.cpp bspline.hpp
	g++ -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o bspline_bind.so bspline_bind.cpp bspline.cpp ${CXXFLAGS} -lboost_python  -lboost_numpy -framework Python

diag.out: diag.o mat.o chkopts
	-${CLINKER} -o $@ diag.o mat.o ${SLEPC_LIB}	

include ${SLEPC_DIR}/conf/slepc_common
