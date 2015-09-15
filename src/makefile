$CXXFLAGS=

pep8:
	pep8 bspline.py

bspline_bind.so: bspline_bind.cpp bspline.cpp bspline.hpp
	g++ -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o bspline_bind.so bspline_bind.cpp bspline.cpp ${CXXFLAGS} -lboost_python  -lboost_numpy -framework Python -Wall

