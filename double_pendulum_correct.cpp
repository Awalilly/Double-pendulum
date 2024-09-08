#include <iostream>  
#include <iomanip>   
#include <cmath>     
#include <fstream>   
#include <stdlib.h>  

using namespace std;

// Function prototypes
void dYdt(double t, double *Y, double *R);  
void RK4Step(double t, double *Y, double h, void (*RHSFunc)(double, double *, double *), int neq);  
double Energy(double theta1, double theta2, double omega1, double omega2);  

double G_ACCg = 9.81;  // Acceleration due to gravity

// Masses and lengths of the pendulums
# define M1 1
# define M2 1
# define L1 1
# define L2 1

// SESSION 1 (Is energy conserved?), SESSION 2 (GIF), SESSION 3 (Chaos), SESSION 4 (Flip)
#define SESSION 1

// PARAMETER == 1 select theta1, PARAMETER == 2 select theta2, PARAMETER == 3 select omega1, PARAMETER == 4 select omega2
#define PARAMETER 4

int main() {

	// Parameters and variables for computation
	double h = 0.;        // Step size 
	double tb = 0.;       // Initial time 
	double t = tb;        // Current time 
	double t_end = 1.e2;  // End time for simulation (100s)
	double E_0, E, relErr; // Energy variables

	// Initial conditions for the double pendulum
	double x1, x2, y1, y2;	// Cartesian coordinates for the masses
	int neq = 4;   	// Number of equations
	int end = 1.e6; // Total number of iterations 
	double Y[neq];	// Initial state vector
    
#if SESSION == 1

	ofstream file; // Output file stream

	file.open("Double_pendulum_energy.txt"); // Open a file to write energy data

	// Declaration and initialization of variables
	double hbeg = 1.e-4; // Initial step size
	double hend = 0.05; // Final step size
	double npoints = 100; // Number of points for step size variation
	double IC[neq]; // Array to hold initial conditions
	
	// Assign initial condition
	Y[0] = 1.6;
	Y[1] = 1.;
	Y[2] = 2.;
	Y[3] = -3.;
	
	// Save the intial conditions
	for(int i = 0; i < neq; i++) IC[i] = Y[i];
	
	// Compute initial energy of the system
	E_0 = Energy(Y[0], Y[1], Y[2], Y[3]);
	E = E_0;

	// Loop for evaluating energy preservation over various step sizes
	for(int j = 1; j <= npoints; j++){

		// Incrementing step size
		h += fabs(hbeg - hend) / npoints;
		t = tb; // Resetting time
        
		// Restoring initial conditions
		for(int k = 0; k < neq; k++) Y[k] = IC[k];
        
		// Performing integration using Runge-Kutta method
		for(int i = 0; i < end; i++){
			RK4Step(t, Y, h, &dYdt, neq);
			t += h; 
		}
        
		t -= h; // Correcting time
        
		// Computing energy and relative error
		E = Energy(Y[0], Y[1], Y[2], Y[3]); 
		relErr = fabs((E_0 - E) / E_0); 
        
		// Writing step size and relative error to file
		file << h << "\t" << relErr << endl;

		// Clear output and print the progression of the process
		system("clear");
		cout << j << endl;
	}   
    
	// Closing the file
	file.close();

#endif
        
#if SESSION == 2

	ofstream file; // Output file stream
	file.open("Motion_small_angles.txt"); // Open a file to write motion data

	// Setting initial conditions for the double pendulum system
	Y[0] = double(5) * M_PI / 180.; // Initial angle of the first pendulum arm in radians
	Y[1] = double(0) * M_PI / 180.; // Initial angle of the second pendulum arm in radians
	Y[2] = 0.; // Initial angular velocity of the first pendulum arm
	Y[3] = 0.; // Initial angular velocity of the second pendulum arm

	h = 0.001; // Selected step size based on the analysis of the previous session
	end = 1.e5; 

	// Iterating over the system dynamics
	for(int i = 0; i < end; i++){
		RK4Step(t, Y, h, &dYdt, neq); 
        
		// Computing positions of pendulum masses
		x1 = L1 * sin(Y[0]); // x-coordinate of the first pendulum mass
		y1 = - L1 * cos(Y[0]); // y-coordinate of the first pendulum mass
		x2 = x1 + L2 * sin(Y[1]); // x-coordinate of the second pendulum mass
		y2 = y1 - L2 * cos(Y[1]); // y-coordinate of the second pendulum mass
        
		// Computing energy and relative error
		E = Energy(Y[0], Y[1], Y[2], Y[3]); 
		relErr = fabs((E_0 - E) / E_0); 
        
		// Writing time, positions, angles, and velocities to file
		file << t << " " << x1 << " " << y1 << " " << x2 << " " << y2 << " "
			 << Y[0] << " " << Y[1] << "\t" << Y[2] << "\t" << Y[3] << endl;

		t += h; // Updating time
	}

	file.close();

#endif


#if SESSION == 3

	// Setting step size
	h = 0.005; // Selected step size based on the analysis of the previous session

	// Array for initial conditions
	double IC[neq];
    
	// Declaration of variables for Lyapunov exponent computation
	double Ly1, Ly2; // Lyapunov exponents
	int n_iter = 100; // Number of iterations on a single trajectory
	double d1, d2, d0_1, d0_2; // Variables for computing distances
	double delta = 1.e-8; // Perturbation
	double x1_a, x2_a, y1_a, y2_a; // Auxiliary variables for pendulum positions
    
	double theta_end = M_PI; // Range for theta
	double omega_end = 2 * M_PI; // Range for omega
	double N_point = 180; // Number of points
	double step; // Step for parameter
	t = 0.; 
	t_end = 1.e3; 
	end = t_end / h; // Total number of iterations
    
	// Setting initial conditions
	Y[0] = 0.;
	Y[1] = 0.;
	Y[2] = 0.;
	Y[3] = 0.;
    
	int param; // The index of the parameter of interest
	double Deviation[] = {delta, delta, 0., 0.}; // Perturbation
    
	// Determining chaos search based on parameter type

#if PARAMETER == 1 // Search chaos for theta1

	ofstream file_chaos;
	file_chaos.open("Theta1_100.txt");
	file_chaos << setiosflags(ios::scientific);
	file_chaos << setiosflags(ios::scientific);
	step = theta_end / N_point;
	param = 0;

#elif PARAMETER == 2 // Search chaos for theta2

	ofstream file_chaos;
	file_chaos.open("Theta2_100.txt");
	file_chaos << setiosflags(ios::scientific);
	file_chaos << setiosflags(ios::scientific);
	step = theta_end / N_point;
	param = 1;

#elif PARAMETER == 3 // Search chaos for Omega1

	ofstream file_chaos;
	file_chaos.open("Omega1_100.txt");
	file_chaos << setiosflags(ios::scientific);
	file_chaos << setiosflags(ios::scientific);
	step = omega_end / N_point;
	param = 2;

#elif PARAMETER == 4 // Search chaos for Omega2

	ofstream file_chaos;
	file_chaos.open("Omega2_100.txt");
	file_chaos << setiosflags(ios::scientific);
	file_chaos << setiosflags(ios::scientific);
	step = omega_end / N_point;
	param = 3;

#endif

	// Looping over parameter values
	for(int i = 0; i <= N_point; i++){ 


		//Set the initial condition
		for(int z = 0; z < neq; z++){
			if(z != param) IC[z] = 0.;
			else IC[z] = step * i;
		}

        
		// Looping for the number of Lyapunov exponents
		for(int j = 0; j < n_iter; j++){

			// Set initial condition with deviation
			for(int k = 0; k < neq; k++) Y[k] = IC[k] + Deviation[k];

			// Compute pendulum positions
			x1 = L1 * sin(IC[0]);
			y1 = - L1 * cos(IC[0]);
			x2 = x1 + L2 * sin(IC[1]);
			y2 = y1 - L2 * cos(IC[1]);

			// Position for deviated coordinates
			x1_a = L1 * sin(Y[0]);
			y1_a = - L1 * cos(Y[0]);
			x2_a = x1_a + L2 * sin(Y[1]);
			y2_a = y1_a - L2 * cos(Y[1]);

			// Initial distance between the two trajectories for both the masses
			d0_1 += sqrt((x1_a - x1) * (x1_a - x1) + (y1_a - y1) * (y1_a - y1));
			d0_2 += sqrt((x2_a - x2) * (x2_a - x2) + (y2_a - y2) * (y2_a - y2));


			// Integrate the deviated system
			for(int k = 0; k < end; k++){
				t += h;
				RK4Step(t, Y, h, &dYdt, neq);
			}

			t = 0.;	// Reset time

			// Compute new pendulum deviated positions
			x1_a = L1 * sin(Y[0]);
			y1_a = - L1 * cos(Y[0]);
			x2_a = x1_a + L2 * sin(Y[1]);
			y2_a = y1_a - L2 * cos(Y[1]);

			// Reset initial condition to the non deviated ones
			for(int k = 0; k < neq; k++) Y[k] = IC[k];


			// Integrate the non deviated the system
			for(int k = 0; k < end; k++){
				t += h;
				RK4Step(t, Y, h, &dYdt, neq);
			}

			// Compute new pendulum non deviated positions
			x1 = L1 * sin(Y[0]);
			y1 = - L1 * cos(Y[0]);
			x2 = x1 + L2 * sin(Y[1]);
			y2 = y1 - L2 * cos(Y[1]);

			t = 0.;	// Reset time

			// Compute the distance at the end of the integration and sum 
			d1 += sqrt((x1_a - x1) * (x1_a - x1) + (y1_a - y1) * (y1_a - y1));
			d2 += sqrt((x2_a - x2) * (x2_a - x2) + (y2_a - y2) * (y2_a - y2));

			// Set the final trajectory values as the new initial conditions
			for(int k = 0; k < neq; k++) IC[k] = Y[k];
		}
        
		// Compute the Lyapunov exponents
		Ly1 = log2(fabs(d1 / d0_1)) / t_end;
		Ly2 = log2(fabs(d2 / d0_2)) / t_end;
        
		// Write to file
		file_chaos << (step * i) * 180 / M_PI << " " << Ly1 << " " << Ly2 << endl;
        
		// Reset distance variables
		d0_1 = 0.;
		d0_2 = 0.;
		d1 = 0.;
		d2 = 0.;

		// Clear output and print the progression of the process
		system("clear");
		cout << i << endl;
	}
    
	// Close file
	file_chaos.close();

#endif


#if SESSION == 4

	// Opening files to store flip times
	ofstream file_flip1, file_flip2;
	file_flip1.open("Flip1.txt");
	file_flip2.open("Flip2.txt");

	// Variables to track if each pendulum has flipped
	bool flipped1 = false;
	bool flipped2 = false;

	// Setting step size and total number of iterations
	h = 0.01;
	end = 1.e5;

	// Initializing velocities
	Y[2] = 0.;
	Y[3] = 0.;

	// Setting parameters for angle computation
	int data_points = 720; // Number of data points
	int deg_lim = data_points / 2; // Range of computation
	double div = 180. / double(deg_lim); // Increment in angle value

	double theta_old1, theta_old2; // Variables to store previous angles

	// Looping through angle combinations
	for(int i = -deg_lim; i <= deg_lim; i++){
		for(int j = -deg_lim; j <= deg_lim; j ++){

			// Setting initial angles
			Y[0] = (double(i) * div) * M_PI / 180.;
			Y[1] = (double(j) * div) * M_PI / 180.;
			Y[2] = 0.;
			Y[3] = 0.;
            
			// Integration loop
			for(int k = 0; k < end; k++){
				RK4Step(t, Y, h, &dYdt, neq);
                
				// Checking if pendulum 1 has flipped
				if((flipped1 == false && t > 1.) 
					&& ((sin(Y[0]) * sin(theta_old1) < 0. 
					&& Y[0] != theta_old1 
					&& cos(Y[0]) < -1. + 1.e-7))){
					
					file_flip1 << (double(i) * div) << "\t" << (double(j) * div) 
										 << "\t" << t << endl; // Writing flip time
								
					flipped1 = true;
				}
                
				if(flipped1 == false) theta_old1 = Y[0]; // Updating previous angle
                
				// Checking if pendulum 2 has flipped
				if((flipped2 == false && t > 1.) 
					&& ((sin(Y[1]) * sin(theta_old2) < 0. 
					&& Y[1] != theta_old2 
					&& cos(Y[1]) < -1. + 1.e-7))){
					
					file_flip2 << (double(i) * div) << "\t" << (double(j) * div) 
										 << "\t" << t << endl; // Writing flip time
					
					flipped2 = true;
				}
                
				// Updating previous angle
				if(flipped2 == false) theta_old2 = Y[1]; 
                
				// Exiting loop if both pendulums have flipped
				if(flipped1 == true && flipped2 == true) break; 
                
				t += h;
			}

			// Writing -1 if pendulum did not flip within the given time
			if(flipped1 == false) {
				file_flip1 << (double(i) * div) << "	" 
									 << (double(j) * div) << "	" << -1 << endl;
			}
			
			if(flipped2 == false){ 
				file_flip2 << (double(i) * div) << "	" 
									 << (double(j) * div) << "	" << -1 << endl;
			}
            
			t = 0.; // Resetting time
			flipped1 = false; // Resetting flip flag for pendulum 1
			flipped2 = false; // Resetting flip flag for pendulum 2
		}
        
		system("clear"); 
		cout << int(i * div) << endl; // Displaying progress
        
		file_flip1 << endl; 
		file_flip2 << endl;
	}
    
	// Closing files
	file_flip1.close();
	file_flip2.close();

#endif

	return 0;
}


void dYdt(double t, double *Y, double *R){
    
	double theta1 = Y[0];
	double theta2 = Y[1];
	double omega1 = Y[2];
	double omega2 = Y[3];
	
	
	double g = G_ACCg;
	double den = 2. * M1 + M2 - M2 * cos(2. * theta1 - 2. * theta2);

	R[0] = omega1;
	R[1] = omega2;
	R[2] = (-g * (2. * M1 + M2) * sin(theta1) - M2 * g * sin(theta1 - 2. * theta2) - 2. * sin(theta1 - theta2) * M2 * (omega2 * omega2 * L2 + omega1 * omega1 * L1 * cos(theta1 - theta2))) / (L1 * den);
	R[3] = (2. * sin(theta1 - theta2) * (omega1 * omega1 * L1 * (M1 + M2) + g * (M1 + M2) * cos(theta1) + omega2 * omega2 * L2 * M2 * cos(theta1 - theta2))) / (L2 * den);
	
}



void RK4Step(double t, double *Y, double h, void (*RHSFunc)(double, double *, double *), int neq){
    
	double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];
    
	RHSFunc(t, Y, k1);
    
	for(int i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * h * k1[i];
    
	RHSFunc(t + 0.5 * h, Y1, k2);
    
	for(int i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * h * k2[i];
    
	RHSFunc(t + 0.5 * h, Y1, k3);
    
	for(int i = 0; i < neq; i++) Y1[i] = Y[i] + h * k3[i];
    
	RHSFunc(t + h, Y1, k4);
    
	for(int i = 0; i < neq; i++) Y[i] += (h / 6.) * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    
    
}

double Energy(double theta1, double theta2, double omega1, double omega2){
	
	return .5 * M1 * L1 * L1 * omega1 * omega1 + .5 * M2 * (L1 * L1 * omega1 * omega1 
			+ L2 * L2 * omega2 * omega2 + 2 * L1 * L2 * omega1 * omega2 * cos(theta1 - theta2)) 
			- (M1 + M2) * G_ACCg * L1 * cos (theta1) - M2 * G_ACCg * L2 * cos (theta2);
}
