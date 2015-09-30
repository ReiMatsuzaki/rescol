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

