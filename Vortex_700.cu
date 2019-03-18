#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

__global__
void NewVortexDistrub(float *V, float *O, float *VN, float *ON, float *S, int N) { // the kernel carries out the assigned formulas on each 
											// thread. i.e. each thread carries the following program.
	__shared__ float Vc[3];			// This kernel calculates angular velocity , Energy, Maximum speed and radius of vortex 
	__shared__ float dVx[3];
	__shared__ float dVy[3];		// Shared variables to give threads access to these variables.
	__shared__ float dVz[3];
	__shared__ float domdt[3];

 int tx = threadIdx.x;			// cuda inbuilit variable threadIdx.x keeps the track of thread id in a block, is assigned to tx
 int bx =  blockIdx.x;			// cuda inbuilt variable blockIdx.x keeps track of blocks assigned during the kernel call
float radiika, dssss_dr;
float ssss;
float t1,t2,t3;

	if (tx<3) {					// Initializes the variable
		if (tx == 0) {
			Vc[tx] = 1.0f;
		}
		else {
			Vc[tx] = 0.0f;
			dVx[tx] = 0.0f;
			dVy[tx] = 0.0f;
			dVz[tx] = 0.0f;
		}
	}

    t1 = ((V[(bx * 3) + 1] - V[(tx * 3) + 1]) * O[(tx * 3) + 2] - (V[(bx * 3) + 2] - V[(tx * 3) + 2]) * O[(tx * 3) + 1]);
    t2 = ((V[(bx * 3) + 2] - V[(tx * 3) + 2]) * O[(tx * 3 + 0)] - (V[(bx * 3) + 0] - V[(tx * 3) + 0]) * O[(tx * 3) + 2]);
    t3 = ((V[(bx * 3) + 0] - V[(tx * 3) + 0]) * O[(tx * 3 + 1)] - (V[(bx * 3) + 1] - V[(tx * 3) + 1]) * O[(tx * 3) + 0]);
	radiika = powf(V[(bx * 3) + 0] - V[(tx * 3) + 0], 2) + pow(V[(bx * 3) + 1] - V[(tx * 3) + 1], 2) + pow(V[(bx * 3) + 2] - V[(tx * 3) + 2], 2);
	ssss = expf((-radiika)*((3.1416f * 2.0f) / (S[tx] * S[tx])));
	dssss_dr = expf(-((3.1416f * 2.0f) / (S[tx] * S[tx])));


	for (int i = 0; i<N; i++) {
		if (tx == i) { 					
			Vc[0] = Vc[0] + ssss* t1;
			Vc[1] = Vc[1] + ssss* t2;
			Vc[2] = Vc[2] + ssss* t3;

			dVx[0] = dssss_dr*(V[(bx * 3) + 0] - V[(tx * 3) + 0]) * t1 + dVx[0];
			dVx[1] = dssss_dr*(V[(bx * 3) + 1] - V[(tx * 3) + 1]) * t1 + O[tx * 3 + 2] * ssss + dVx[1];
			dVx[2] = dssss_dr*(V[(bx * 3) + 2] - V[(tx * 3) + 2]) * t1 - O[tx * 3 + 1] * ssss + dVx[2];

			dVy[0] = dssss_dr*(V[(bx * 3) + 0] - V[(tx * 3) + 0]) * t2 - O[tx * 3 + 2] * ssss + dVy[0];
			dVy[1] = dssss_dr*(V[(bx * 3) + 1] - V[(tx * 3) + 1]) * t2 + dVy[1];
			dVy[2] = dssss_dr*(V[(bx * 3) + 2] - V[(tx * 3) + 2]) * t2 + O[tx * 3 + 0] * ssss + dVy[2];

			dVz[0] = dssss_dr*(V[(bx * 3) + 0] - V[(tx * 3) + 0]) * t3 + O[tx * 3 + 1] * ssss + dVz[0];
			dVz[1] = dssss_dr*(V[(bx * 3) + 1] - V[(tx * 3) + 1]) * t3 - O[tx * 3 + 0] * ssss + dVz[1];
			dVz[2] = dssss_dr*(V[(bx * 3) + 2] - V[(tx * 3) + 2]) * t3 + dVz[2];
		}
		__syncthreads(); // barrier set to get all the threads to this point, before any further calculations
	}

	if (tx<3) {
		VN[(bx * 3) + tx] = V[(bx * 3) + tx] + 0.01f * Vc[tx];  // translates to global thread ids performing the tasks.
		domdt[tx] = dVx[0] * O[bx * 3 + 0] + dVx[1] * O[bx * 3 + 1] + dVx[2] * O[bx * 3 + 2];
		ON[bx * 3 + tx] = O[bx * 3 + tx] + domdt[tx] * 0.01f;
	}
}


__global__
 void AssignNewToOld(float* V, float* O, float* VN, float* ON, float* S, float* A, float* vx, float* E, float* S_m, float* Ss, int N) {
	//*A = 0.0f;  			// This kernel calculates the magnitude of angular velocity , Energy, Maximum speed and radius of vortex.
	//*vx = 0.0f;
	float Amagnit_old;
	float Amagnit_new;
	int tx = threadIdx.x;
	__shared__ float v[3];
	__shared__ float radiika;
	__shared__ float t1;
	__shared__ float Om22P;
	__shared__ float ssss;

	V[(tx * 3) + 0] = VN[(tx * 3) + 0];
	V[(tx * 3) + 1] = VN[(tx * 3) + 1];
	V[(tx * 3) + 2] = VN[(tx * 3) + 2];

	Amagnit_old = sqrtf(powf(O[tx * 3 + 0], 2) + powf(O[tx * 3 + 1], 2) + powf(O[tx * 3 + 2], 2));
	O[tx * 3 + 0] = ON[tx * 3 + 0];
	O[tx * 3 + 1] = ON[tx * 3 + 1];
	O[tx * 3 + 2] = ON[tx * 3 + 2];
	Amagnit_new = sqrtf(powf(O[tx * 3 + 0], 2) + powf(O[tx * 3 + 1], 2) + powf(O[tx * 3 + 2], 2));
	S[tx] = S[tx] * sqrtf(Amagnit_old / Amagnit_new);


	for (int i = 0; i<N; i++) {
		if (tx == i) {
			if (Amagnit_new >= (*A)) {
				*A = Amagnit_new;
				*E = (pow(Amagnit_new, 2)* pow(S[tx], 5));
				*S_m = Amagnit_new * S[tx];
				*Ss = S[tx];
			}
		}
		__syncthreads();
	}

	for(int i=0;i<N;i++){
		if(tx == i){
			v[0] = 0.5 - V[tx * 3 + 0];
			v[1] = 0.5 - V[tx * 3 + 1];
			v[2] = 0.5 - V[tx * 3 + 2];
			radiika = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
			t1 = v[1] * O[tx * 3 + 2] - v[2] * O[tx * 3 + 1];
			Om22P = (3.1416*2.0) / (S[tx] * S[tx]);
			ssss = expf(-radiika*Om22P);
			*vx = *vx + (ssss*t1);
		}
		__syncthreads();
	}
}

float randomgenerator()    		//Returns a random number between 0 and 1
{
	float x;
   	x=rand()/float(RAND_MAX);	    //Rand_MAX value depends on library, but is always less than 32767 on any standard library usage
   	return x;
}

float randomgenerator(float a, float b)  //Returns a random number between the passed parameters
{
	 float t;
     t=(b-a)*randomgenerator() + a;   //The values a and b (a>b) are the range in which the numbers are required.
	 return t;
}


int main(int argc, char ** argv) {

	clock_t start, end;
	start = clock();		// For time calculation
	int Ntime = 10000;        // Number of time steps
	float Delta_t = 0.01;  // Time step
	float Radius = 0.1;   // Radius of vorton (the same for all)
	int number = 700;                           // number of vortons
	int Ncout;
	float Replace;
	float StatisticalMoments[4];

	for (int i = 0; i<4; i++) {
		StatisticalMoments[i] = 0.0;
	}
	//--------------Define Host and device arrays------------------------//
	
	static float h_Vortex[2100];			// Arrays are defined for host and devices
	static float h_Omega_v[2100];
	static float h_VortexN[2100];
	static float h_Omega_vN[2100];
	static float h_Sigma[700];
	float h_Amagni;
	float h_Energy;
	float h_Speed_max;
	float h_Sigmas;
	float h_Vx;
	int size1 = sizeof(float);

	float *d_Vortex;
	float *d_Omega_v;
	float *d_VortexN;
	float *d_Omega_vN;
	float *d_Sigma;
	float *d_Amagni;
	float *d_Energy;
	float *d_Speed_max;
	float *d_Sigmas;
	float *d_Vx;

	cudaMalloc((void **)&d_Vortex, number * 3 * size1);			// Size of variable is allocated on GPU's dram
	cudaMalloc((void **)&d_Omega_v, number * 3 * size1);
	cudaMalloc((void **)&d_VortexN, number * 3 * size1);
	cudaMalloc((void **)&d_Omega_vN, number * 3 * size1);
	cudaMalloc((void **)&d_Sigma, number*size1);
	cudaMalloc((void **)&d_Amagni, size1);
	cudaMalloc((void **)&d_Energy, size1);
	cudaMalloc((void **)&d_Speed_max, size1);
	cudaMalloc((void **)&d_Sigmas, size1);
	cudaMalloc((void **)&d_Vx, size1);

	/*---------------------------------------Creating files------------------------------------------------*/

	//ofstream outfile1;
	//ofstream outfile2;
	// outfile1.open("./outputfiles/velocities_cuda.dat");
	// outfile2.open("./outputfiles/MaxValue_cuda.dat");

	/* ----------------------------------------Intilization of Arrays---------------------------------------*/

	for (int i = 0; i<number * 3; i++) { 			// Host arrays are initialized with a random function, giving them 	
		srand(i);									// values in the range of (0, 1)
		h_Vortex[i] = randomgenerator(1,0);
		srand(i);
		h_Omega_v[i] =randomgenerator(1,0)- 0.5;	
	}

	for (int i = 0; i<number; i++) {
		srand(i);
		h_Sigma[i] = Radius;
	}

	//----------------------------Transfer Data of host to Device -------------------//

	cudaMemcpy(d_Vortex, h_Vortex, number * 3 * size1, cudaMemcpyHostToDevice);	
	cudaMemcpy(d_Omega_v, h_Omega_v, number * 3 * size1, cudaMemcpyHostToDevice);	
	cudaMemcpy(d_Sigma, h_Sigma, number*size1, cudaMemcpyHostToDevice);				

	//--------------------------------------------------------------------------------------------------------//

	for (int itime = 0; itime<Ntime; itime++) {

		NewVortexDistrub << <number, number >> >(d_Vortex, d_Omega_v, d_VortexN, d_Omega_vN, d_Sigma, number); // kernel call with 
													// 700 blocks each with 700 threads.
		cudaDeviceSynchronize();					// All threads required to hit a barrier before any further calculations.
		cudaMemcpy(h_Omega_vN, d_Omega_vN, number * 3 * size1, cudaMemcpyDeviceToHost); // Data transfer back to host array
		cudaMemcpy(h_VortexN, d_VortexN, number * 3 * size1, cudaMemcpyDeviceToHost);	// Data transfer back to host array

		/*--------------------------mapping to the cube back--------------------------------------*/
		Ncout = 0;
		for (int ivorton = 0; ivorton< number; ivorton++) {
			Replace = 0.0f;
			for (int i = 0; i <3; i++) {
				if (h_VortexN[ivorton * 3 + i]< 0.0f) { Replace = 1.0f; }
				if (h_VortexN[ivorton * 3 + i] > 1.0f) { Replace = 1.0f; }
			}
			if (Replace == 1.0) {
				Ncout = Ncout + 1;
				srand(ivorton);// intialization of co-ordinate of vortrons
				h_VortexN[ivorton * 3 + 0] = randomgenerator(1,0);
				h_VortexN[ivorton * 3 + 1] = randomgenerator(1,0);
				h_VortexN[ivorton * 3 + 2] = randomgenerator(1,0);
				// intilization of strength of vortrons
				h_Omega_vN[ivorton * 3 + 0] = randomgenerator(1,0) - 0.5;
				h_Omega_vN[ivorton * 3 + 1] = randomgenerator(1,0) - 0.5;
				h_Omega_vN[ivorton * 3 + 2] = randomgenerator(1,0) - 0.5;
				h_Sigma[ivorton] = Radius;   // intilization of radius of vortrons
			}
		}

		cudaMemcpy(d_Omega_vN, h_Omega_vN, number * 3 * size1, cudaMemcpyHostToDevice);
		cudaMemcpy(d_VortexN, h_VortexN, number * 3 * size1, cudaMemcpyHostToDevice);

		/*------------------------- mapping to the cube back---------------------------------------*/
		//the old parameters became new ones
		h_Amagni = 0.0f;
		h_Vx = 0.0f;
		cudaMemcpy(d_Amagni,&h_Amagni,size1,cudaMemcpyHostToDevice);
		cudaMemcpy(d_Vx,&h_Vx,size1,cudaMemcpyHostToDevice);
		AssignNewToOld <<<1, number >>>(d_Vortex, d_Omega_v, d_VortexN, d_Omega_vN, d_Sigma, d_Amagni, d_Vx, // 2nd kernel call for 
d_Energy, d_Speed_max, d_Sigmas, number);			// calculating magnitude of the arrays passed as parameters
		cudaDeviceSynchronize();

		cudaMemcpy(&h_Amagni, d_Amagni, size1, cudaMemcpyDeviceToHost);			// Data Transfers made accordingly
		cudaMemcpy(&h_Energy, d_Energy, size1, cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_Speed_max, d_Speed_max, size1, cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_Sigmas, d_Sigmas, size1, cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_Vx, d_Vx, size1, cudaMemcpyDeviceToHost);


		printf("%15.3f%15.3f%15.3e%15.3e%15.3e%15.3f", itime*Delta_t, h_Amagni, h_Energy, h_Speed_max, h_Sigmas, h_Vx);
		printf("\n");
		for (int ier = 0; ier<4; ier++) {
			StatisticalMoments[ier] = StatisticalMoments[ier] + pow(h_Vx, ier);
		}
	}

	end = clock();
	printf("Time taken for GPU code for 700X700 threads is %10.3f sec" , (float) (end-start)/CLOCKS_PER_SEC ); //compute time
	cudaFree(d_Vortex);					// memory freed after the calculations
	cudaFree(d_Omega_v);
	cudaFree(d_VortexN);
	cudaFree(d_Omega_vN);
	cudaFree(d_Sigma);
	cudaFree(d_Amagni);
	cudaFree(d_Vx);
	cudaFree(d_Energy);
	cudaFree(d_Speed_max);
	cudaFree(d_Sigmas);
	
	return 0;
}


