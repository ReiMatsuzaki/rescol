WORK=work
LMAX=4
FEM="-fem_type bss -bss_order 5 -bps_num_zs 51 -bps_type exp -bps_zmax 50.0"
EPS="-eps_nev 5 -eps_max_it 1000 -eps_type krylovschur"
BOND="-bondlength0 0.0 -num_bondlength 21 -d_bondlength 0.1"

mkdir -p $WORK
rm ${WORK}/*

python ../../../script/calc_y2mat.py -t ${WORK} -l1 ${LMAX} > y2mat.dat

../../../src_c/he_guess.out -in_dir ${WORK} -out_dir ${WORK} ${FEM} ${EPS} \
			    -z 1.5  -guess_type calc | tee he_guess.dat && \
../../../src_c/h2mole.out -eri direct -guess_type read  -qmax 32 \
			  -in_dir ${WORK} -out_dir ${WORK} ${FEM} ${EPS} \
			  ${BOND} -malloc_log | tee h2mole.dat 
	
 
