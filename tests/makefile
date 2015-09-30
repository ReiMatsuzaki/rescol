export PYTHONPATH=../src

test_bind: test_bind.py
	python test_bind.py

test_bspline: test_bspline.py
	pep8 test_bspline.py
	pep8 ../src/bspline.py
	python test_bspline.py

test_hmole: test_hmole.py
	pep8 test_hmole.py
	pep8 ../src/hmole.py
	python test_hmole.py

test_fdiff: test_fdiff.py
	pep8 test_fdiff.py
	pep8 ../src/fdiff.py
	python test_fdiff.py

test_utils: test_utils.py
	pep8 test_utils.py
	pep8 ../src/utils.py
	python test_utils.py

test_angmoment: test_angmoment.py
	python --version
	pep8 test_angmoment.py
	pep8 ../src/angmoment.py
	python test_angmoment.py

test_all: test_bind test_bspline test_hmole test_utils test_angmoment
