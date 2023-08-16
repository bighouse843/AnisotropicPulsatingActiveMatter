/*

  This code simulates the dynamics of N active deforming particles:
  
  \dot r_i = - \mu_r \sum_{jvi} \nabla U(r_i-r_j, \phi_i, \phi_j) + \sqrt{2\mu_r T} \xi_i
	\dot\phi_i = \omega + \mu_\phi \sum_{jvi} [ \epsilon sin(\phi_i-\phi_j) - \partial_{\phi_i} U(r_i-r_j, \phi_i, \phi_j)] + \sqrt{2\mu_\phi T} \eta_i

	where

	range of interaction    = radius * [1 + lambda * sin(\phi)] / (1 + lambda)
	< \xi_i(t) \xi_j(0) >   = \delta_{ij} \delta(t)
	< \eta_i(t) \eta_j(0) > = \delta_{ij} \delta(t)


	Lx, Ly          = system size
	rho             = density
  T               = temperature
  mu_r            = mobility of position
  mu_phi          = mobility of size variable
  mu_theta        = mobility of the angle variable
  omega           = driving of size variable
  dt              = time step
  totaltime       = duration of the run
  initial_time    = waiting time before record
  interval_record = time interval between successive records
  interval_snap   = time interval between successive snapshots
  radius          = range of interaction
  amplitude_r     = strength of interaction
  amplitude_phi   = strength of synchronization
  amplitude_theta   = strength of alignment
	lambda          = interaction parameter
	dist            = radius of density averaging: 2*dist+radius < min(Lx,Ly)
	dx              = binning of density correlations
	drho            = binning of density
	dphi            = binning of size

*/
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include <vector>
#include "mtwister.c"
#include <random>
#include <deque>

#define PI 3.14159265358979323846

double U = 0.0;

using namespace std;
// Structure Particule: contains coordinates and box id of each particle
typedef struct {
  double qx; // Position
  double qy;
  double phi; // Phase (governing minor axes length)
  double theta; // Orientation
} particle;

// Structure Forces: contains x and y projection of interaction forces with neighboring particles
typedef struct {
  double fx;
  double fy;
  double torque_phi;
  double torque_theta;
} force;

// class for a site object, containing a list of particles on it
class site{

private:
  vector<int> particle_ids;
  
public:

  vector<int> get_particle_ids(){
    return particle_ids;
  }

  void remove_particle_id(int id){
    for (int k = 0; k < this->particle_ids.size(); ++k){
      if(particle_ids[k]==id){
	particle_ids.erase(particle_ids.begin() + k);
      }
    }
  }

  void copy_particle_ids_to_deque(deque<int> &d){
    for(auto id : this->particle_ids){
      d.push_back(id);
    }
  }

  
  void add_particle_id(int id){
    this->particle_ids.push_back(id);
  }

};

// site = linked chained of all particle in a given box
// elem: itself + pointer to the next particles
typedef struct elem elem;

struct elem
{
  int num; // number of a particle in the box
  struct elem *next; // pointer onto the next element. next->num is the number of another particle in the box.
};
	
//typedef elem* site; // Site is a pointer on a box, i.e. on a particle and on a pointer to another particle.


// Compute the sign of an argument
int sign(double x){
  return (x > 0) ? 1 : (-1);
}
int sign_int(int x){
  return (x > 0) ? 1 : (-1);
}

double floor0(double x){
  return (x>=0) ? floor(x) : ceil(x);
}


// Choose next and previous elements with periodic boundary condition
int NEXT(int k, double L){
  return (k == (int) L-1) ? 0 : (k+1);
}

int PREV(int k, double L){
  return (k == 0) ? ((int) L-1) : (k-1);
}

int IDX(int k, double L, int n){
  return (k == (int) L+n) ? 0 : (k+n);
}


/* ---------------------------
    DECLARATION OF FUNCTIONS
   ------------------------------*/

// Store data
void Record(FILE* outputdata, vector<particle> &Particles, int N); // Particle coordinates


// Start from particles randomly deposited in the system
void InitialConditions(int N, int Nx, int Ny, double Lx, double Ly, double radius, double a, vector<vector<site>> &sites, vector<particle> &Particles, MTRand &ran);

// Compute forces between neighboring particles
void add_force(int k1, particle p1, int k2, particle p2, double p2x, double p2y, double radius, double lambda, double amplitude_r, double amplitude_phi, double amplitude_theta, double Ly, double Lx, vector<force> &Force);

// Implement the dynamics
double MoveParticles(int N, int Nx, int Ny, double Lx, double Ly, vector<double> &Noise_r, vector<double> &Noise_phi, vector<double> &Noise_theta, double omega, double radius, double lambda, double amplitude_r, double amplitude_phi, double amplitude_theta, double mu_r, double mu_phi, double mu_theta, vector<vector<site>> &sites,vector<particle> &Particles, double dt, double Time, vector<force> &Force);



/* ---------------------------
    END DECLARATION OF FUNCTIONS
   ------------------------------*/

/* ---------------------------
    MAIN
   ------------------------------*/

int main(int argc, char *argv[]){ 
  /* 
     argc is the number of argument of the command line, and argv is a list
     of string containing the arguments. argv[0] is the name of the
     executable.
  */
  
  char params[] = "Input: output Lx Ly rho T mu_r mu_phi mu_theta omega dt totaltime initial_time interval_record interval_snap seed radius amplitude_r amplitude_phi amplitude_theta lambda dist drho dx dphi dtheta";

  // Check that the number of inputs is correct
  if(argc!=26){
    printf("%s\n",params);
    exit(1);
  }
  // ./out test 4 4 0.125 0.001 1 1 1 0.01 0.1 100 1 1 1 1234 1 1 1 1 0.3 0.001 0.1 0.1 0.1 0.1
  //g++ -o out Param_order.cpp


  // DEFINITION OF VARIABLES  
  int i, j, k, idx; // Iterators
  int N; // Number of particles
  int Nx, Ny; // Number of boxes of length radius in Lx, Ly
  int idx_frame = 0; // Index of the frame
  double Lx, Ly; // System size
  double rho; // Particle density
  double T; // Temperature
  double mu_r, mu_phi, mu_theta; // Mobilities
  double omega; // Driving of size
  double amplitude_r, amplitude_phi, amplitude_theta; // Strength of interactions
	double drho, dx, dphi, dtheta, rho_max, dist; // Binning of histogram
	long   N_corr, N_dens, N_phi, N_theta; // Maximum value of histograms
	int    idx_max; // Size of density averaging in box length
	double mu_avg_phi,mu2_avg_phi,mu_avg_theta,mu2_avg_theta;
	double Size_phi, Size2_phi, Current_phi, dummy_phi; // Order parameter and current of size
	double Size_theta, Size2_theta, Current_theta, dummy_theta; // Order parameter and current of size
  double dt; // Time step
  double Time; // Current time
  double totaltime; // Total duration of the run (after initialisation)
  double initial_time; // Time for equilibration
  double interval_record; // Time between two recordings
  double interval_snap; // Time between two snapshots
  double nb_record;  // Number of time steps of the next record
  double nb_snap;  // Number of time steps of the next snapshot
  long   counter_stat, aux; // Counter for statistics
  long   seed; // Seed of random number generator
  double radius, lambda; // Particle radius
  double a;
  time_t clock; // Time to measure duration of the simulation
  
  FILE *output; // File where parameters are stored
  FILE *output_snap; // File where snapshots are stored
  FILE *output_corr; // File where correlation is stored
  FILE *output_dens; // File where density is stored
  FILE *output_phi; // File where size is stored
  FILE *output_theta; // File where orientation is stored
  FILE *output_order_phi; // File where size order parameter is stored
  FILE *output_order2_phi;
  FILE *output_order_theta; // File where size order parameter is stored
  FILE *output_order2_theta;
  FILE *output_order_avg;
  FILE *output_order2_avg;
  char name[200]; // Name of output file containing the parameters
  char name_snap[215]; // Name of output file containing snapshots
  char name_corr[205]; // Name of output file containing correlation
  char name_dens[205]; // Name of output file containing density
  char name_phi[205]; // Name of output file containing size
  char name_theta[210]; // Qui il 205 dava errore, aggiornato a 210
  char name_order_phi[210]; // Name of output file containing order parameter
  char name_order2_phi[215];
  char name_order_theta[215]; // Name of output file containing order parameter
  char name_order2_theta[215];


  // READING THE INPUT
  i = 1;
  sprintf(name, "%s", argv[i]);
  sprintf(name_corr, "%s_corr", name);
  sprintf(name_dens, "%s_dens", name);
  sprintf(name_phi, "%s_phi", name);
  sprintf(name_theta, "%s_theta", name);
  sprintf(name_order_phi, "%s_order_phi", name);
  sprintf(name_order2_phi, "%s_order2_phi", name);
  sprintf(name_order_theta, "%s_order_theta", name);
  sprintf(name_order2_theta, "%s_order2_theta", name);
  output      = fopen(name, "w");
  output_order_phi = fopen(name_order_phi, "w");
  output_order2_phi = fopen(name_order2_phi, "w");
  output_order_theta = fopen(name_order_theta, "w");
  output_order2_theta = fopen(name_order2_theta, "w");

  i++;
  Lx              = strtod(argv[i], NULL); i++;
  Ly              = strtod(argv[i], NULL); i++;
  rho             = strtod(argv[i], NULL); i++;
  T               = strtod(argv[i], NULL); i++;
  mu_r            = strtod(argv[i], NULL); i++;
  mu_phi          = strtod(argv[i], NULL); i++;
  mu_theta          = strtod(argv[i], NULL); i++;
  omega           = strtod(argv[i], NULL); i++;
  dt              = strtod(argv[i], NULL); i++;
  totaltime       = strtod(argv[i], NULL); i++;
  initial_time    = strtod(argv[i], NULL); i++;
  interval_record = strtod(argv[i], NULL); i++;
  interval_snap   = strtod(argv[i], NULL); i++;
  seed            = strtol(argv[i], NULL, 10); i++;
  radius          = strtod(argv[i], NULL); i++;
  amplitude_r     = strtod(argv[i], NULL); i++;
  amplitude_phi   = strtod(argv[i], NULL); i++;
  amplitude_theta   = strtod(argv[i], NULL); i++;
  lambda          = strtod(argv[i], NULL); i++;
  dist            = strtod(argv[i], NULL); i++;
  drho            = strtod(argv[i], NULL); i++;
  dx              = strtod(argv[i], NULL); i++;
  dphi            = strtod(argv[i], NULL); i++;
  dtheta            = strtod(argv[i], NULL); i++;

  // Print the parameters in output  
  fprintf(output, "%lg\n", Lx);
  fprintf(output, "%lg\n", Ly);
  fprintf(output, "%lg\n", rho);
  fprintf(output, "%lg\n", T);
  fprintf(output, "%lg\n", mu_r);
  fprintf(output, "%lg\n", mu_phi);
  fprintf(output, "%lg\n", mu_theta);
  fprintf(output, "%lg\n", omega);
  fprintf(output, "%lg\n", dt);
  fprintf(output, "%lg\n", totaltime);
  fprintf(output, "%lg\n", initial_time);
  fprintf(output, "%lg\n", interval_record);
  fprintf(output, "%lg\n", interval_snap);
  fprintf(output, "%lg\n", radius);
  fprintf(output, "%lg\n", amplitude_r);
  fprintf(output, "%lg\n", amplitude_phi);
  fprintf(output, "%lg\n", amplitude_theta);
  fprintf(output, "%lg\n", lambda);
  fprintf(output, "%lg\n", dist);
  fprintf(output, "%lg\n", drho);
  fprintf(output, "%lg\n", dx);
  fprintf(output, "%lg\n", dphi);
  fprintf(output, "%lg\n", dtheta);
  fflush(output);


  // INITIALISATION OF VARIABLES
 // init_genrand64(seed);
  MTRand ran = seedRand(seed);
  // Number of particles given rho, Lx and Ly
  a = sqrt(2/(sqrt(3)*rho));
  N =  static_cast<int>(2*round(Lx/a)*round(Ly/(sqrt(3)*a)));
  

  // Number of boxes of length radius in Lx, Ly
  //printf("2radius %.3f\n", 2*radius);
  //printf("Lx/2radius %.3f\n", Lx/(2*radius));
  //printf("floor0(Lx/2radius) %.3f\n", floor0(Lx/(2*radius)));
  // printf("static_cast<int>floor0(Lx/2radius) %.d\n", static_cast<int>(floor0(Lx/(2*radius))));
  Nx = static_cast<int>(floor(Lx/(2*radius*exp(lambda))));
  Ny = static_cast<int>(floor(Ly/(2*radius*exp(lambda))));  

  if(Nx<3){
    printf("Nx < 3 !");
    exit(1);
  }
  if(Ny<3){
    printf("Ny < 3 !");
    exit(1);
  }

  // Start the clock time
  clock = time(NULL);

  // Initially all linked chains are NULL pointers
   vector<vector<site>> sites(Nx,vector<site>(Ny));

  // Initialize arrays
  vector<double> Noise_r(2*N,0); // noise array (thermal)
  vector<double> Noise_phi(N,0); // noise array (active)
  vector<double> Noise_theta(N,0); // noise array (thermal)
  particle p; // generic starter particle
  vector<particle> Particles(N,p); // particles
  force f; // generic starter force
  vector<force> Force(N,f); // forces on the particles

  // Initial conditions
  InitialConditions(N, Nx, Ny, Lx, Ly, radius, a, sites, Particles,ran);
  printf("Random deposition succeeded\n");
  Time      = -initial_time;
  Size_phi         = 0;
  Current_phi      = 0;
  Size_theta         = 0;
  Current_theta      = 0;
  nb_record    = 0;
  nb_snap      = 0;
  Size2_phi		= 0;
  Size2_theta		= 0;
  mu_avg_phi	= 0;
  mu2_avg_phi	= 0;
  mu_avg_theta	= 0;
  mu2_avg_theta	= 0;
	counter_stat = 0;
	idx          = 0;
  output_order_phi = fopen(name_order_phi, "w");
  output_order_theta = fopen(name_order_theta, "w");

  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<double> d{0.,1.};
   // Run dynamics until the time iterator reaches the final time
    
    
  while(Time<totaltime){

    // Sample realizations of the noise 
    for (k=0;k<N;k++){
    	Noise_r[2*k]   = sqrt(2*mu_r*T)*(d(gen));
			Noise_r[2*k+1] = sqrt(2*mu_r*T)*(d(gen));
			Noise_phi[k]   = sqrt(2*mu_phi*T)*(d(gen));
      Noise_theta[k]   = sqrt(2*mu_theta*T)*(d(gen));
		}

    // Move the particles according to the dynamics
	 Time = MoveParticles(N, Nx, Ny, Lx, Ly, Noise_r, Noise_phi, Noise_theta, omega, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, mu_r, mu_phi, mu_theta, sites, Particles, dt, Time, Force);

    // After thermalization
    if(Time>=nb_snap){
	 		// Print simulation progress
			//printf("%s: %.4lg\n", name, 1e2*Time/totaltime);

			// Allow record of statistics
			idx = 1;

			// Increase snap counter
      
				nb_snap += interval_snap;
			

/*
				 if(idx==1){

                  // Compute order parameter and current of size
                  dummy_phi = 0;
                  Size_phi = 0;
                  Size2_phi = 0.;
                  dummy_theta = 0;
                  Size_theta = 0;
                  Size2_theta = 0.;
                  // Current=0;
                    for(i=0;i<N;i++){
                        //	 Current += omega + mu_phi*Force[i].torque;
                          for(j=i;j<N;j++)   {   
                            dummy_phi += cos(2*Particles[i].phi-2*Particles[j].phi);
                            dummy_theta += cos(2*Particles[i].theta-2*Particles[j].theta);
                          }
                        }
                    Size_phi += sqrt(2*dummy_phi);
                    Size2_phi += 2*dummy_phi;
                    Size_theta += sqrt(2*dummy_theta);
                    Size2_theta += 2*dummy_theta;

                    // Increment of statistics counter
                    counter_stat++;
               		 	}

        //fprintf(output_order_phi, "%.15lg\t%.15lg\n",Time+initial_time, U);
				fprintf(output_order_phi, "%.15lg\t%.15lg\n",Time+initial_time, Size_phi/N);
				fflush(output_order_phi);
				fprintf(output_order2_phi, "%.15lg\t%.15lg\n",Time+initial_time, Size2_phi/(N*N));
        fflush(output_order2_phi);
				mu_avg_phi += Size_phi/N;
				mu2_avg_phi += Size2_phi/(N*N);
				fprintf(output_order_theta, "%.15lg\t%.15lg\n",Time+initial_time, Size_theta/N);
				fflush(output_order_theta);
				fprintf(output_order2_theta, "%.15lg\t%.15lg\n",Time+initial_time, Size2_theta/(N*N));
        fflush(output_order2_theta);
				mu_avg_theta += Size_theta/N;
				mu2_avg_theta += Size2_theta/(N*N);
*/

		}

			if(Time>=nb_record){
        idx_frame++;
        printf("%d\n", idx_frame);
        sprintf(name_snap, "%s_snap_%d", name, idx_frame);
        output_snap  = fopen(name_snap, "w");
				// Record snapshot of system
				Record(output_snap, Particles, N);

				// Increase snap counter
        nb_record += interval_record;
				
			}
  }
	
	
  // Return the duration of the simulation
  printf("Simulation duration: %ld seconds\n", (long) time(NULL)-clock);
  fprintf(output, "%ld\n", (long) time(NULL)-clock);
  fprintf(output, "%ld\n", counter_stat);
  fprintf(output, "%f\n", mu_avg_phi);
  fprintf(output, "%f\n", mu2_avg_phi);
  fprintf(output, "%f\n", mu_avg_theta);
  fprintf(output, "%f\n", mu2_avg_theta);
  fclose(output_order_phi);
  fclose(output_order_theta);
  fclose(output);
  return 0;
}

/* ---------------------------
    END MAIN
   ------------------------------*/

/* ---------------------------
    DYNAMICS
   ------------------------------*/

double MoveParticles(int N, int Nx, int Ny, double Lx, double Ly, vector<double> &Noise_r, vector<double> &Noise_phi, vector<double> &Noise_theta, double omega, double radius, double lambda, double amplitude_r, double amplitude_phi, double amplitude_theta, double mu_r, double mu_phi, double mu_theta, vector<vector<site>> &sites, vector<particle> &Particles, double dt, double Time, vector<force> &Force){
  
  int k1, k2, k, i, j, i_new, j_new, next_i, prev_i, prev_j;
  double dqx, dqy, dphi, dtheta; // Increment of coordinates
	double test, max_amp=0; // Adaptative time stepping
	double dt_g, sdt_g; // Time step
  double temp_x, temp_y;

  particle p1, p2;
  site site1; // Iterator on the particles of the site (i,j)
  site site2; // Iterator on the particles of the neighboring sites of the site (i,j)
   for(i=0;i<N;i++){
    Force[i].fx=0.0;
    Force[i].fy=0.0;
    Force[i].torque_phi=0.0;
    Force[i].torque_theta=0.0;
  } 
 U=0;

  deque<int> temp_ids;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      // create a copy of particle ids to utilize newton's third law (correctly!) on same site                    
      sites[i][j].copy_particle_ids_to_deque(temp_ids);
      next_i = NEXT(i,Nx);
      prev_i = PREV(i,Nx);
      prev_j = PREV(j,Ny);
      
        		for(auto k1: sites[i][j].get_particle_ids()){ // Iteraction on particles in site (i,j)
				
				// pop the first element off temp_ids (using deque container ensures operation is order 1 (not n))
        		temp_ids.pop_front();

				p1 = Particles[k1]; // Structure of particle 1
			for(auto k2: temp_ids ){ // Loop over the particle ids
				  
				  p2 = Particles[k2]; // Get its structure
					// Compute forces from particles on the SAME site
					add_force(k1, p1, k2, p2, Particles[k2].qx, Particles[k2].qy, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, Ly, Lx, Force); 
					//site2 = site2->next;
				}
			// Site on the right
			for(auto k2: sites[next_i][j].get_particle_ids()){ // Loop over the particle ids
				 p2 = Particles[k2]; // Get its structure
          temp_x = Particles[k2].qx;
          temp_y = Particles[k2].qy;
          
            //printf("i: %d \n", i);
            //printf("j: %d \n", j);
            //printf("next_i: %d \n", next_i);
          if (next_i == 0) {
            temp_x = Particles[k2].qx + Lx;
          }

          // Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, temp_x, temp_y, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, Ly, Lx, Force); 
				}
	
				
				
			// Site on the bottom left (NOTA CHE I E J SONO ORIENTATE COME X E Y)
			for(auto k2: sites[prev_i][prev_j].get_particle_ids()){ // Loop over the particle ids	 
				  p2 = Particles[k2]; // Get its structure
          temp_x = Particles[k2].qx;
          temp_y = Particles[k2].qy;

          if (prev_j == (int) Ny-1 ) {
            temp_y = Particles[k2].qy - Ly;
          }
          if (prev_i == (int) Nx-1) {
            temp_x = Particles[k2].qx - Lx;
            //printf("bottom left if temp y: %.3f \n", temp_y);
          }

          // Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, temp_x, temp_y, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, Ly, Lx, Force); 
				}
	
				
		// Site on the bottom 
		for(auto k2: sites[i][prev_j].get_particle_ids()){ // Loop over the particle ids		 
				  p2 = Particles[k2]; // Get its structure
          temp_x = Particles[k2].qx;
          temp_y = Particles[k2].qy;

          if (prev_j == (int) Ny-1) {
            temp_y = Particles[k2].qy - Ly;
            //printf("bottom if temp y: %.3f \n", temp_y);
          }


          // Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, temp_x, temp_y, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, Ly, Lx, Force); 
				}
	
				
		// Site on the bottom right
		for(auto k2: sites[next_i][prev_j].get_particle_ids()){ // Loop over the particle ids				  
				  p2 = Particles[k2]; // Get its structure
					
          temp_x = Particles[k2].qx;
          temp_y = Particles[k2].qy;

          
          if (prev_j == (int) Ny-1 ) {
            temp_y = Particles[k2].qy - Ly;
          }
          if (next_i == 0) {
            temp_x = Particles[k2].qx + Lx;
            //printf("top left if temp y: %.3f \n", temp_y);
          }

          // Compute forces from neighboring particles
					add_force(k1, p1, k2, p2, temp_x, temp_y, radius, lambda, amplitude_r, amplitude_phi, amplitude_theta, Ly, Lx, Force); 
				}
      }
    }
  }

	// Compute maximum of interaction amplitude
	for(k=0;k<N;k++){
		test = mu_r*fabs(Force[k].fx)/(2*radius);
		if(test>max_amp) max_amp = test;
		test = mu_r*fabs(Force[k].fy)/(2*radius);
		if(test>max_amp) max_amp = test;
		test = (omega + mu_phi*Force[k].torque_phi)/(2*PI);
		if(test>max_amp) max_amp = test;
    test = (mu_theta*Force[k].torque_theta)/(2*PI);
		if(test>max_amp) max_amp = test;
	}



 	// Adaptative time stepping
	if(max_amp*dt<1e-1){
		dt_g = dt;
		} else{
	dt_g = 1e-1/max_amp;
	}
	sdt_g    = sqrt(dt_g);
	Time += dt_g;	
  // Update coordinates of all particles
  for(k=0;k<N;k++){ // Iterations over the N particles
    
    // Extract position and label of the box containing the particle
    i = static_cast<int>(floor(Particles[k].qx*Nx/Lx)); 
    j = static_cast<int>(floor(Particles[k].qy*Ny/Ly));

	  // Dynamics
		dqx  = dt_g*(mu_r*Force[k].fx) + sdt_g*Noise_r[2*k]; //TOGLI ADDENDO 1*K
	  dqy  = dt_g*(mu_r*Force[k].fy) + sdt_g*Noise_r[2*k+1]; //TOGLI ADDENDO 1*K
		dphi = dt_g*(omega + mu_phi*Force[k].torque_phi) + sdt_g*Noise_phi[k];
    dtheta = dt_g*(mu_theta*Force[k].torque_theta) + sdt_g*Noise_theta[k];


    // New coordinates after one iteration
    Particles[k].qx  += dqx;
    Particles[k].qy  += dqy;
    Particles[k].phi += dphi;
    Particles[k].theta += dtheta;


    // Pay attention to the boundary conditions
    if(Particles[k].qx>Lx || Particles[k].qx<0)  Particles[k].qx -= Lx*floor(Particles[k].qx/Lx);
		if(Particles[k].qy>Ly || Particles[k].qy<0)  Particles[k].qy -= Ly*floor(Particles[k].qy/Ly);

    // New label of the box containing the particle
    i_new = static_cast<int>(floor(Particles[k].qx*Nx/Lx)); 
    j_new = static_cast<int>(floor(Particles[k].qy*Ny/Ly));
    //printf("k: %d \n",k);
    //printf("i_new: %d \n",i_new);
    //printf("j_new: %d \n",j_new);

    //printf("k: %d \n\n",k);

    // Update de Sites
    if(i!=i_new || j!=j_new){
     sites[i_new][j_new].add_particle_id(k);
      sites[i][j].remove_particle_id(k);
    }
  }
  return Time;
}

/* ---------------------------
    END DYNAMICS
   ------------------------------*/


/* ---------------------------
    FORCES
   ------------------------------*/

// Compute the force between particles and add it to the force array
void add_force(int k1, particle p1, int k2, particle p2, double p2x, double p2y, double radius, double lambda, double amplitude_r, double amplitude_phi, double amplitude_theta, double Ly, double Lx, vector<force> &Force){

  double sigma, sigma6, rx, ry, r, r6, u1x, u1y, u2x, u2y, r12x, r12y, l1, l2, l1_2, l2_2; // Interparticle distance
  double amp_rep, amp_phi, amp_theta, amp_align; // Interaction strengths
  double sigma0, d1, d2, chi2, alpha2, tau, tau_, sigma0_2, d1_2, d2_2, radius_2, dr_dx, dr_dy, dpsi_dx, dpsi_dy, l1_2_m_d1_2, l2_2_p_d1_2, l1_2_p_d2_2, l2_2_m_d2_2;
  double r12u1, r12u2, u1u2, r12u1_2, r12u2_2, u1u2_2, sigma3, dU_dr, dU_dpsi, dU_dx, dU_dy, dU_dsigma, dsigma_dpsi, psi, dd1_dphi, dd2_dphi2, dsigma_dsigma1, dsigma1_dr12u1, dsigma1_du1u2, dsigma1_dr12u2;
  double tau_r12u1_2, btau_r12u2_2, chi2u1u2, r12u1r12u2, chi2_1, r12pu1, r12pu2, tau_r12u1, btau_r12u2, dsigma_dtheta, sigma3_sigma0_2, dsigma_dd1, dsigma_dd2, dl2_dphi2, dl1_dphi, d1_2_p_d2_2;
  double chi2_1_2, u1pu2, r12u1p, dU_dtheta, dsigma_dsigma0, lsin_i_1, lcos_i, lsin_i_1_2, lsin_j_1, lsin_j_1_2, l_1, l_1_2, sinj_1, lsinj_l_2, lsinj_l_2_l_1_2_2, sinj;
  double dsigma0_dphi, lcos_i_lsin_i_1, dtau_dphi, lsinj_1_2_l_1_2, dbtau_dphi, dchi2_dphi, lsini_1_2_l_1_2, lsini_1_2_l_1_2_2, dsigma_dtau, dsigma_dbtau, dsigma_dchi2, dU_dphi;
  double dtau_dphi2, dbtau_dphi2, lcos_j_lsin_j_1, lcos_j, lsinj_1_2_l_1_2_2, lsini_l_2, sini, sini_1, dchi2_dphi2, dsigma0_dphi2, k0_1_p_l, k0_1_m_l, dU_dphi2, u1u2p, r12u2p, dsigma_dtheta2, dU_dtheta2;
  
  // Interparticle distance
  rx  =  p1.qx - p2x;
  ry  =  p1.qy - p2y;
  //printf("temp p2x: %.3f \n", p2x);
  //printf("temp p2y: %.3f \n", p2y);
  //printf("p2x: %.3f \n", p2.qx);
  //printf("p2y: %.3f \n", p2.qy);
  //printf("rx: %.3f \n", rx);
  //printf("ry: %.3f \n", ry);

	r = sqrt(rx*rx + ry*ry);
  r12x = rx/r;
  r12y = ry/r;
  u1x = cos(p1.theta);
  u1y = sin(p1.theta);
  u2x = cos(p2.theta);
  u2y = sin(p2.theta);
  k0_1_p_l = radius/(1+lambda); //k0/(1+lambda)
  k0_1_m_l = radius*(1-lambda); //k0*(1-lambda)
  d1 = radius *exp(lambda*sin(p1.phi));
  d2 = radius *exp(lambda*sin(p2.phi));
  d1_2 = d1*d1;
  d2_2 = d2*d2;
  l1 = radius *exp(-lambda*sin(p1.phi));
  l2 = radius *exp(-lambda*sin(p2.phi));
  l1_2 = l1*l1;
  l2_2 = l2*l2;
  sigma0 = sqrt(2)*sqrt(d1_2+d2_2);
  l1_2_m_d1_2 = l1_2 - d1_2;
  l2_2_p_d1_2 = l2_2 + d1_2;
  l1_2_p_d2_2 = l1_2 + d2_2;
  l2_2_m_d2_2 = l2_2 - d2_2;
  chi2 = (l1_2_m_d1_2*l2_2_m_d2_2)/(l2_2_p_d1_2*l1_2_p_d2_2);
	tau = l1_2_m_d1_2/l1_2_p_d2_2;
  tau_ = l2_2_m_d2_2/l2_2_p_d1_2;
  dr_dx = r12x; //segni corretti, questo e successivo
  dr_dy = r12y;
  r12u1 = r12x*u1x+r12y*u1y;
  r12u2 = r12x*u2x+r12y*u2y;
  u1u2 = u1x*u2x+u1y*u2y;
  r12u1_2 = r12u1*r12u1;
  r12u2_2 = r12u2*r12u2;
  u1u2_2 = u1u2*u1u2;
  tau_r12u1_2 = tau*r12u1_2; //A
  btau_r12u2_2 = tau_*r12u2_2; //B
  chi2u1u2 = chi2*u1u2; //C
  r12u1r12u2 = r12u1*r12u2; //D
  chi2_1 = 1-chi2*u1u2_2; //1-chi^2(u1u2)^2
  sigma = sigma0/sqrt(1-(tau_r12u1_2+btau_r12u2_2-2*chi2u1u2*r12u1r12u2)/chi2_1);


	// Add contribution to the force
	if(r<sigma){
    //printf("r: %.3f \n", r);
    //printf("sigma: %.3f \n", sigma);
    //printf("Interazione \n");
		// Interaction strengths
  //printf("r: %.3f \n", r);
   // printf("sigma: %.3f \n", sigma);
    dpsi_dx = -dr_dy/r;
    dpsi_dy = +dr_dx/r;
    sigma3 = sigma*sigma*sigma;
    sigma6 = sigma3*sigma3;
    r6 = r*r*r*r*r*r;
    sigma0_2 = sigma0*sigma0;
    r12pu1 = -r12y*u1x+r12x*u1y;
    r12pu2 = -r12y*u2x+r12x*u2y;
    tau_r12u1 = tau*r12u1;
    btau_r12u2 = tau_*r12u2;
    sigma3_sigma0_2 = sigma3/sigma0_2;


    //Amplitude_r è epsilon0
    //Position component of the force
    dU_dr = 12*sigma6*amplitude_r/r6/r*(1-sigma6/r6); //verifica il segno
    dU_dsigma = -dU_dr*r/sigma;
    
    dsigma_dsigma1 = sigma3_sigma0_2/2;

    dsigma1_dr12u1  = 2*(tau_r12u1-chi2u1u2*r12u2)/chi2_1;
    dsigma1_dr12u2 = 2*(btau_r12u2-chi2u1u2*r12u1)/chi2_1;

    dsigma_dpsi = dsigma_dsigma1*(dsigma1_dr12u1*r12pu1+dsigma1_dr12u2*r12pu2);
    dU_dpsi = dU_dsigma*dsigma_dpsi;
		
    dU_dx = +dU_dr*dr_dx+dU_dpsi*dpsi_dx;
    dU_dy = +dU_dr*dr_dy+dU_dpsi*dpsi_dy;

    U += amplitude_r*(sigma6*sigma6/(r6*r6)-2*sigma6/r6+1);
   //printf("dU: %.3f \n",U.Utot);

    //Orientation component of the force acting on the first particle
    chi2_1_2 = chi2_1*chi2_1; //(1-chi^2(u1u2)^2)^2
    u1pu2 = -u1y*u2x+u1x*u2y;
    r12u1p = -r12x*u1y+u1x*r12y;
    
    dsigma1_du1u2 = 2*chi2*(u1u2*(btau_r12u2_2 + tau_r12u1_2) - r12u1r12u2*(2-chi2_1))/chi2_1_2;

    dsigma_dtheta = dsigma_dsigma1*(dsigma1_dr12u1*r12u1p+dsigma1_du1u2*u1pu2);
    
    dU_dtheta = dU_dsigma*dsigma_dtheta;


    //Phase component of the force acting on the first particle
    lcos_i = lambda*cos(p1.phi);
    dd1_dphi = lcos_i * d1;
    dl1_dphi = - lcos_i*l1;
    d1_2_p_d2_2 = d1_2 + d2_2;
 

    dsigma_dsigma0 = sigma/sigma0;
    dsigma0_dphi = 2/sigma0*d1*dd1_dphi; 

    dsigma_dtau = r12u1_2;//già raccolto sigma3_sigma0_2/(2*chi2_1)
    dtau_dphi = -2/l1_2_p_d2_2*d1*dd1_dphi + 2*dl1_dphi*l1*d1_2_p_d2_2/(l1_2_p_d2_2*l1_2_p_d2_2); 

    dsigma_dbtau = r12u2_2;//già raccolto sigma3_sigma0_2/(2*chi2_1)
    dbtau_dphi = -2*tau_/l2_2_p_d1_2*d1*dd1_dphi; 

    dsigma_dchi2 = (tau_r12u1_2*u1u2_2 + btau_r12u2_2*u1u2_2-2*r12u1r12u2*u1u2)/chi2_1;//già raccolto sigma3_sigma0_2/(2*chi2_1)
    dchi2_dphi = tau*dbtau_dphi + tau_*dtau_dphi; 

    dsigma_dd1 = dsigma_dsigma0*dsigma0_dphi+ sigma3_sigma0_2/(2*chi2_1)*(dsigma_dtau*dtau_dphi+dsigma_dbtau*dbtau_dphi+dsigma_dchi2*dchi2_dphi);

    dU_dphi = dU_dsigma*dsigma_dd1;

//printf("lcos_i: %.3f \n", lcos_i); 
//printf("dd1_dphi: %.3f \n", dd1_dphi ); 

//printf("dU_dphi: %.3f \n", dU_dphi); 




    //Phase component of the force acting on the second particle
    lcos_j = lambda*cos(p2.phi);
    dd2_dphi2 = lcos_j * d2;
    dl2_dphi2 = - lcos_j*l2;

    dsigma0_dphi2 = 2/sigma0*d2*dd2_dphi2; 

    dtau_dphi2 = -2/l1_2_p_d2_2*d2*dd2_dphi2*tau; //-tau*2*d2/(l1^2 + d2^2)

    dbtau_dphi2 = -2/l2_2_p_d1_2*d2*dd2_dphi2 + 2*dl2_dphi2*l2*d1_2_p_d2_2/(l2_2_p_d1_2*l2_2_p_d1_2);

    dchi2_dphi2 = tau*dbtau_dphi2 + tau_*dtau_dphi2;

    dsigma_dd2 = dsigma_dsigma0*dsigma0_dphi2+ sigma3_sigma0_2/(2*chi2_1)*(dsigma_dtau*dtau_dphi2+dsigma_dbtau*dbtau_dphi2+dsigma_dchi2*dchi2_dphi2);

    dU_dphi2 = dU_dsigma*dsigma_dd2;

    //Orientation component of the force acting on the second particle
    r12u2p = -r12x*u2y+u2x*r12y;
    u1u2p = - u1pu2;

    dsigma_dtheta2 = dsigma_dsigma1*(dsigma1_dr12u2*r12u2p+dsigma1_du1u2*u1u2p);
 
    dU_dtheta2 = dU_dsigma*dsigma_dtheta2;

		amp_align = amplitude_phi*sin(p2.phi-p1.phi);		
    //amp_theta

		// Interaction on p1
		Force[k1].fx     -= dU_dx;
		Force[k1].fy     -= dU_dy;
		Force[k1].torque_phi += amp_align*cos(p2.phi-p1.phi) - dU_dphi;
    Force[k1].torque_theta -= dU_dtheta; //verifica bene questo segno
    //printf("Force theta particle %d: ", k1);
    //printf("Force: %.3f \n", Force[k1].torque_phi ); 
		
		// Interaction on p2 COME SI SCRIVE PRECISAMENTE?
		Force[k2].fx     += dU_dx;
		Force[k2].fy     += dU_dy;
		Force[k2].torque_phi -= amp_align*cos(p1.phi-p2.phi) + dU_dphi2; 
    Force[k2].torque_theta -= dU_dtheta2; //verifica bene questo segno
    //printf("Force theta particle %d: ", k2);
    //printf("%.3f \n", Force[k2].torque_theta );
	}
}

/* ---------------------------
    END FORCES
   ------------------------------*/

/* ---------------------------
    INITIAL CONDITIONS
   ------------------------------*/

// Deposition of N particles uniformly in the system
void InitialConditions(int N, int Nx, int Ny, double Lx, double Ly, double radius, double a, vector<vector<site>> &sites, vector<particle> &Particles,MTRand &ran){


  int k = 0, k1 = 0, k2 = 0; // Iterator on particles
  int i, j; // Indices of the box
  double X, Y; // Sampled coordinates
  double shift = a/3; 

  // Iterations over the deposited particles
  while(k<N){
    
    if ((k1*a + shift) > Lx) {
      k1 = 0;
      k2++;
    }
    if ((k1*a + (k2%2)*a/2) > Lx ) {
      k1 = 0;
      k2++;
    }
    X = k1*a + (k2%2)*a/2 + shift;
    Y = k2*a*sqrt(3)/2;
    k1++;
    if(X>Lx || X<0)  X -= Lx*floor(X/Lx);
		if(Y>Ly || Y<0)  Y -= Ly*floor(Y/Ly);


    //if(k==0){
    //  X = Lx/2;
    //  Y = Ly/2;
    //}
    //else{
    //  X = 0;
    //  Y = Ly/2;
    //}

    // Extract label (i,j) of corresponding box
    i = static_cast<int>(floor0(X*Nx/Lx));
    j = static_cast<int>(floor0(Y*Ny/Ly));

    //printf("i: %d \n",i);
    //printf("j: %d \n",j);
  
    // Store spatial coordinates
    Particles[k].qx = X;
    Particles[k].qy = Y;

  
    // Initial size variables
    Particles[k].phi = -PI/2; //originale
    Particles[k].theta = PI/6;



    //if(k==0){
    // Particles[k].phi = 0;
    //  Particles[k].theta = 0;
    //}
    //else{
    //  Particles[k].phi = 0;
    //  Particles[k].theta = 0;     
    //}    
 
    // Add particle in the system
    sites[i][j].add_particle_id(k);
    k++;
  }
}

/* ---------------------------
    END INITIAL CONDITIONS
   ------------------------------*/


/* ---------------------------
    PRINT
   ------------------------------*/

// Record particle coordinates
void Record(FILE* outputdata, vector<particle> &Particles, int N){

  int i;
  
  for(i=0;i<N;i++){
    fprintf(outputdata, "%lg\t%lg\t%lg\t%lg\n", Particles[i].qx, Particles[i].qy, Particles[i].phi, Particles[i].theta);
  }

  fflush(outputdata);
}

/* ---------------------------
    END PRINT
   ------------------------------*/


