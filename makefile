main: setup
	mpicxx -I. -I${CMF}/include -I${PTL}/include -O3 main.cc -o pbc -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

setup:
	mkdir -p output
	mkdir -p series
	mkdir -p checkpoint

clean:
	rm -f pbc
	rm -r output
	rm -r series
