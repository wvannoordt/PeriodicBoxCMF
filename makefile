FLAGS := -Ofast
export DIM=3
export DOLATEXOUTPUT=0
export OPTLEVEL=3
export PARALLEL=1
export CUDA_ENABLE=0
#FLAGS := -O3
main: setup cmf
	mpicxx --std=c++11 -I. -I${CMF}/include -I${PTL}/include ${FLAGS} main.cc -o pbc -L${CMF}/lib -lcmf -L${PTL}/lib -lPTL

cmf:
	make -C ${CMF} -f makefile

setup:
	mkdir -p output
	mkdir -p series
	mkdir -p checkpoint

clean:
	rm -f pbc
	rm -r output
	rm -r series
