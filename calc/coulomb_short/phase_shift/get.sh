echo "" > ene_phase.dat
grep energy res_8_500.dat | cut -f2 -d: > enes.dat
grep phase res_8_500.dat | cut -f2 -d: > phase.dat
paste enes.dat phase.dat > ene_phase.dat
rm enes.dat phase.dat

#for order in 4 6 8
#do
#    for num in 200 300 500
#    do
#	res="res_"${order}_${num}.dat
#	echo $res
#	grep phase $res | cut -f2 -d: > phase.dat
#	paste ene_phase.dat phase.dat > tmp.dat
#	cp tmp.dat ene_phase.dat
#	rm tmp.dat
#    done
#done

