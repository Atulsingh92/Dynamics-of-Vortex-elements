#include<iostream>
#include<stdio.h>
#include<random>
#include<math.h>
#include<time.h>
#include<fstream>

using namespace std;

double randomgenerator()   		//Returns a random number between 0 and 1
{   
	double x;
   	x=rand()/double(RAND_MAX);	//Rand_MAX value depends on library, but is always less than 32767 on any standard library usage
   	return x;
}

double randomgenerator(double a, double b)//Returns a random number between the passed parameters
{
	 double t;
     t=(b-a)*randomgenerator() + a; //The values a and b (a>b) are the range in which the numbers are required.
	 return t;
}

int main(){
	clock_t start, end;
	start = clock();
	int Ntime=10000	;                               // Number of time steps
	double Delta_t=0.01	;				           // Time step
	double Radius=0.1	;							   // Radius of vorton (the same for all)
	int Number=700	;							   // number of vortons
    double V_mean = 1.0;
	double Amagni;
	double Amagnit_new;
	double Amagnit_old;
	double Energy;
	int Ncout;
	  double Vxc;
	  double Vyc;
	  double Vzc;
	      double dvxdxmov;
	      double dvxdymov;
	      double dvxdzmov;
	      double dvydxmov;
	      double dvydymov;
	      double dvydzmov;
	      double dvzdxmov;
	      double dvzdymov;
	      double dvzdzmov;
		  double vxx;
		  double vyy;
		  double vzz;
		  double radiika;
		  double t1;
		  double t2;
		  double t3;
		  double Om22P;
		  double ssss;
		  double dssss_dr;
		  double domxdt;
		  double domydt;
		  double domzdt;
		  double Replace;
		  double Speed_max;
		  double Sigmas;
		  double Vx;

	double StatisticalMoments[4];
		for(int i=0;i<4;i++){
			StatisticalMoments[i] = 0.0;
		}

    static double Vortex[2100];
    static double Omega_v[2100];
    static double VortexN[2100];
    static double Omega_vN[2100];
    static double Sigma[700];

/* ----------------------------------------Intilization of Arrays---------------------------------------*/
for(int i=0;i<Number*3;i++){
        
			// intialization of co-ordinate of vortrons
			Vortex[i] = randomgenerator(1,0);

			// intilization of strength of vortrons
			Omega_v[i] = randomgenerator(1,0)-0.5;
}

for(int i =0;i<Number;i++){
    Sigma[i] = Radius;  	// intilization of radius of vortrons
}

//---------------------------------------Creating files------------------------------------------------*/

  ofstream outfile1;
  ofstream outfile2;
  outfile1.open("/home/atul/Desktop/cuda_project/files/velocities700.dat"); //Change these paths to the folder in which the source file
  outfile2.open("/home/atul/Desktop/cuda_project/files/maxvalues700.dat");  // is kept. Result shall be obtained in the same folders.

//--------------------------------------------------------------------------------------------------------//

for(int itime=0;itime<Ntime;itime++){

		for (int ivorton = 0; ivorton<Number; ivorton++) { //parallelize this loop
										//calculation of the velocity and tensor S_ij induced at vorton with number "ivorton"
		Vxc = 1.0;
		Vyc = 0.0;
		Vzc = 0.0;
		dvxdxmov = 0.0;
		dvxdymov = 0.0;
		dvxdzmov = 0.0;
		dvydxmov = 0.0;
		dvydymov = 0.0;
		dvydzmov = 0.0;
		dvzdxmov = 0.0;
		dvzdymov = 0.0;
		dvzdzmov = 0.0;

		for (int induced = 0; induced<Number; induced++) { //parallelize this loop
			vxx = Vortex[ivorton * 3 + 0] - Vortex[induced * 3 + 0];
			vyy = Vortex[ivorton * 3 + 1] - Vortex[induced * 3 + 1];
			vzz = Vortex[ivorton * 3 + 2] - Vortex[induced * 3 + 2];

			radiika = vxx*vxx + vyy*vyy + vzz*vzz;
			t1 = vyy*Omega_v[induced * 3 + 2] - vzz*Omega_v[induced * 3 + 1];
			t2 = vzz*Omega_v[induced * 3 + 0] - vxx*Omega_v[induced * 3 + 2];
			t3 = vxx*Omega_v[induced * 3 + 1] - vyy*Omega_v[induced * 3 + 0];
			Om22P = (3.1416 * 2.0) / (Sigma[induced] * Sigma[induced]);
			ssss = exp(-radiika*Om22P);

			//Here are the velocities
			Vxc = Vxc + ssss*t1;
			Vyc = Vyc + ssss*t2;
			Vzc = Vzc + ssss*t3;

		//Here are the strain rate tensor S_ij

			dssss_dr = (-Om22P)*ssss;

			dvxdxmov = dssss_dr*vxx*t1 + dvxdxmov;
			dvxdymov = dssss_dr*vyy*t1 + Omega_v[induced * 3 + 2] * ssss + dvxdymov;
			dvxdzmov = dssss_dr*vzz*t1 - Omega_v[induced * 3 + 1] * ssss + dvxdzmov;

			dvydxmov = dssss_dr*vxx*t2 - Omega_v[induced * 3 + 2] * ssss + dvydxmov;
			dvydymov = dssss_dr*vyy*t2 + dvydymov;
			dvydzmov = dssss_dr*vzz*t2 + Omega_v[induced * 3 + 0] * ssss + dvydzmov;

			dvzdxmov = dssss_dr*vxx*t3 + Omega_v[induced * 3 + 1] * ssss + dvzdxmov;
			dvzdymov = dssss_dr*vyy*t3 - Omega_v[induced * 3 + 0] * ssss + dvzdymov;
			dvzdzmov = dssss_dr*vzz*t3 + dvzdzmov;
		}

		//new coordinates of the vorton "ivorton"	calculated using explicit Euler method
		VortexN[ivorton * 3 + 0] = Vortex[ivorton * 3 + 0] + Delta_t*Vxc;
		VortexN[ivorton * 3 + 1] = Vortex[ivorton * 3 + 1] + Delta_t*Vyc;
		VortexN[ivorton * 3 + 2] = Vortex[ivorton * 3 + 2] + Delta_t*Vzc;

		//new strengths of the vorton "ivorton"  calculated using explicit Euler method
		domxdt = dvxdxmov* Omega_v[ivorton * 3 + 0] + dvxdymov * Omega_v[ivorton * 3 + 1] + dvxdzmov* Omega_v[ivorton * 3 + 2];

		domydt = dvydxmov* Omega_v[ivorton * 3 + 0] + dvydymov* Omega_v[ivorton * 3 + 1] + dvydzmov* Omega_v[ivorton * 3 + 2];

		domzdt = dvzdxmov* Omega_v[ivorton * 3 + 0] + dvzdymov* Omega_v[ivorton * 3 + 1] + dvzdzmov* Omega_v[ivorton * 3 + 2];

		Omega_vN[ivorton * 3 + 0] = Omega_v[ivorton * 3 + 0] + domxdt*Delta_t;
		Omega_vN[ivorton * 3 + 1] = Omega_v[ivorton * 3 + 1] + domydt*Delta_t;
		Omega_vN[ivorton * 3 + 2] = Omega_v[ivorton * 3 + 2] + domzdt*Delta_t;
	}



/*--------------------------mapping to the cube back--------------------------------------*/
Ncout =0;
    for (int ivorton = 0; ivorton< Number; ivorton++) {
			Replace = 0.0;
			for (int i = 0; i <3; i++) {
				if (VortexN[ivorton * 3 + i]< 0.0) { Replace = 1.0; }
				if (VortexN[ivorton * 3 + i] > 1.0) { Replace = 1.0; }
			}
			if (Replace == 1.0) {
				Ncout = Ncout + 1;

  
			// intialization of co-ordinate of vortrons
			VortexN[ivorton * 3 + 0] = randomgenerator(1,0);
			VortexN[ivorton * 3 + 1] = randomgenerator(1,0);
			VortexN[ivorton * 3 + 2] = randomgenerator(1,0);
 
			// intilization of strength of vortrons
			Omega_vN[ivorton * 3 + 0] = randomgenerator(1,0)-0.5;
			Omega_vN[ivorton * 3 + 1] = randomgenerator(1,0)-0.5;
			Omega_vN[ivorton * 3 + 2] = randomgenerator(1,0)-0.5;
			Sigma[ivorton] = Radius; 								// intilization of radius of vortrons
        }
    }
/*------------------------- mapping to the cube back---------------------------------------*/
    Amagni =0.0;
	//the old parameters became new ones
	for (int ivorton = 0; ivorton<Number; ivorton++) {

			Vortex[ivorton * 3 + 0] = VortexN[ivorton * 3 + 0];
			Vortex[ivorton * 3 + 1] = VortexN[ivorton * 3 + 1];
			Vortex[ivorton * 3 + 2] = VortexN[ivorton * 3 + 2];

			Amagnit_old = sqrt(pow(Omega_v[ivorton * 3 + 0], 2) + pow(Omega_v[ivorton * 3 + 1], 2) + pow(Omega_v[ivorton * 3 + 2], 2));
			Omega_v[ivorton * 3 + 0] = Omega_vN[ivorton * 3 + 0];
			Omega_v[ivorton * 3 + 1] = Omega_vN[ivorton * 3 + 1];
			Omega_v[ivorton * 3 + 2] = Omega_vN[ivorton * 3 + 2];

			Amagnit_new = sqrt(pow(Omega_v[ivorton * 3 + 0], 2) + pow(Omega_v[ivorton * 3 + 1], 2) + pow(Omega_v[ivorton * 3 + 2], 2));
			Sigma[ivorton] = Sigma[ivorton] * sqrt(Amagnit_old / Amagnit_new);

			if (Amagnit_new >= Amagni) {
				Amagni = Amagnit_new;
				Energy = pow(Amagnit_new, 2)*pow(Sigma[ivorton], 5);
				Speed_max = Amagnit_new * Sigma[ivorton];
				Sigmas = Sigma[ivorton];
			}
		}

		Vx = 0.0;
		for (int induced = 0; induced<Number; induced++) {
			vxx = 0.5 - Vortex[induced * 3 + 0];
			vyy = 0.5 - Vortex[induced * 3 + 1];
			vzz = 0.5 - Vortex[induced * 3 + 2];
			radiika = vxx*vxx + vyy*vyy + vzz*vzz;
			t1 = vyy*Omega_v[induced * 3 + 2] - vzz*Omega_v[induced * 3 + 1];
			Om22P = (3.1416 * 2.0) / (Sigma[induced] * Sigma[induced]);
			ssss = exp(-radiika*Om22P);
			Vx = Vx + (ssss*t1);
		}

		outfile1 << itime*Delta_t << "       " << Amagni << "       " << Energy << "        " << Speed_max << "          " << Sigmas << endl;
		outfile2 << itime*Delta_t << "        " <<Vx << endl;


    for(int ier=0;ier<4;ier++){
      StatisticalMoments[ier]=StatisticalMoments[ier]+pow(Vx,ier);
	}
}
	end = clock();
	cout <<" Time taken for CPU code for 700X700 threads is   " << (float) (end-start)/CLOCKS_PER_SEC << "s" << endl;
	outfile1.close();
    outfile2.close();
	 return 0;
  }


