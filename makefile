run: main
	mpirun -np 3 ./program

main:
	mpicxx -I. -I${CMF}/include -I${PTL}/include -DCMF_DIM=3 -O3 main.cc -o program -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

clean:
	rm -f program
