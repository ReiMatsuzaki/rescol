res=res_8_101.dat
echo "" > ene_phase.dat
grep energy ${res} | cut -f2 -d: > enes.dat
grep phase ${res} | cut -f2 -d: > phase.dat
grep "abs_amp " ${res} | cut -f2 -d: > abs_amp.dat
paste enes.dat phase.dat > ene_phase.dat
paste enes.dat abs_amp.dat > ene_abs_amp.dat
rm enes.dat phase.dat abs_amp.dat

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

