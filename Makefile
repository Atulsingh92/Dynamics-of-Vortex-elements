#target: dependency
#		command
all: queue
	
queue:	Vortex_700.out
		qsub vortex12new1.job

Vortex_700.out: Vortex_700.cu
				nvcc -arch sm_13 Vortex_700.cu -o Vortex_700.out
