rm -r common/
mkdir common/
cd common
python ../../../../script/calc_y2mat.py -l0 0 -l1 2 -l1guess 0 -l2guess 0 | tee calc_y2mat.out
cd ..

mkdir -p 0
../../../src_c/he_guess.out -L1 0 -L2 0 \
			    -z 1.00 \
			    -in_dir common \
			    -out_dir 0 \
			    -bss_knots_type exp \
			    -bss_rmax 30.0 \
			    -bss_order 2\
			    -bss_knots_num 21 | tee he_guess.out

cp common/* 0 
mkdir -p 1
../../../src_c/h2mole.out -in_dir 0 \
			  -out_dir 1 \
			  -bss_order 2 \
			  -bss_knots_num 21 \
			  -bss_knots_type exp \
			  -bss_rmax 30.0 \
			  -eri direct \
			  -qmax 10 \
			  -bond_length 0.0 \
			  -guess_type vec \
			  -eps_target -4.2 \
			  -eps_nev 1 \

cp common/* 1			  
cp 0/eig0.vec.dat 1/
mv 1/eig0.vec.dat 1/guess.vec.dat
../../../src_c/h2mole.out -in_dir 1 \
			  -out_dir 1 \
			  -bss_order 2 \
			  -bss_knots_num 21 \
			  -bss_knots_type exp \
			  -bss_rmax 30.0 \
			  -eri direct \
			  -qmax 10 \
			  -bond_length 1.0 \
			  -guess_type vec \
			  -eps_target -4.2 \
			  -eps_nev 1 \
