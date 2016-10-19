grep R0 res_r0* | cat | cut -d: -f3  > r0.dat
grep phase res_r0* | cat | cut -d: -f3 > phase.dat
paste r0.dat phase.dat > r0_phase.dat
#rm r0.dat
#rm phase.dat
