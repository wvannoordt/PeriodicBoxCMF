run: main
	mpirun -np 8 ./program

main: setup
	mpicxx -I. -I${CMF}/include -I${PTL}/include -DCMF_DIM=3 -O3 main.cc -o program -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

setup:
	mkdir -p output

clean:
	rm -f program
	rm -r output