echo "# order num phase_shift" > res.dat
for order in 4 6 8
do
    for num in 200 300 500
    do
	res="res_"${order}_${num}.dat
	echo "${order} ${num}" > tmp1.dat
	grep phase $res | cut -f2 -d: > tmp2.dat
	paste tmp1.dat tmp2.dat >> res.dat
    done
done
rm tmp1.dat tmp2.dat

