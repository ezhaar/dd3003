default:
	@echo -e "please try: make all...\n"

usage:
	@echo -e "\nUSAGE: \n\t./jacobi_serial size iterations\n\t./gauss_serial size iterations\n\t./jacobi_parallel size numThreads iterations\n\t./gauss_parallel size numThreads iterations\n\t./multigrid_serial size iterations\n\t./multigrid_gauss_serial size iterations\n\t./multigrid_parallel size numThreads iterations\n\t./multigrid_gauss_parallel size numThreads iterations\n"

all: jacobi multigrid gauss multigrid_gauss usage

jacobi: jacobi_serial jacobi_parallel
	
gauss: gauss_serial gauss_parallel

multigrid: multigrid_serial multigrid_parallel 

multigrid_gauss: multigrid_gauss_serial multigrid_gauss_parallel

jacobi_serial: jacobi_serial.c
	gcc jacobi_serial.c time.c -lpthread -w -o jacobi_serial

jacobi_parallel: jacobi_parallel.c
	gcc jacobi_parallel.c time.c -lpthread -w -o jacobi_parallel

gauss_serial: gauss_serial.c    
	gcc gauss_serial.c time.c -lpthread -w -o gauss_serial

gauss_parallel: gauss_parallel.c
	gcc gauss_parallel.c time.c -lpthread -w -o gauss_parallel

multigrid_serial: multigrid_serial.c
	gcc multigrid_serial.c time.c -lpthread -w -o multigrid_serial

multigrid_parallel: multigrid_parallel.c
	gcc multigrid_parallel.c time.c -lpthread -w -o multigrid_parallel

multigrid_gauss_serial: multigrid_gauss_serial.c
	gcc multigrid_gauss_serial.c time.c -lpthread -w -o multigrid_gauss_serial

multigrid_gauss_parallel: multigrid_gauss_parallel.c
	gcc multigrid_gauss_parallel.c time.c -lpthread -w -o multigrid_gauss_parallel


clean: 
	rm jacobi_serial jacobi_parallel multigrid_serial multigrid_parallel multigrid_gauss_serial multigrid_gauss_parallel gauss_serial gauss_parallel
