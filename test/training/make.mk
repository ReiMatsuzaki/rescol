train_viewer.out : train_viewer.o
check_train_viewer: train_viewer.out
	./train_viewer.out -my_view ascii:tmp/myfile.dat:ascii_info
	cat tmp/myfile.dat
	./train_viewer.out -my_view binary:tmp/myfile.dat:ascii_info

train_pf.out: train_pf.o
check_train_pf: train_pf.out
	./train_pf.out -malloc_dump
