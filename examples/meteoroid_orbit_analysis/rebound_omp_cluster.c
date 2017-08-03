/**
 * OpenMP example.
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator and direct summation.
 * Shared memory parallelization using OpenMP 
 * is enabled in the Makefile.
 *
 * Note that you need a compiler which supports 
 * OpenMP to run this example. By default, the 
 * OSX compilers from Apple do currently not
 * support OpenMP. You can install the GNU 
 * C compilers easily with homebrew. Look at the
 * Makefile of this example to see how you can setup
 * the parameters to compile REBOUND on both OSX
 * and Linux.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
//#include <string.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"

// Planet Radii (km)
r_mercury = 2439;               
r_venus = 6052;                 
r_earth = 6378.1363;            
r_mars = 3397.3;                
r_jupiter = 71492;              
r_saturn = 60268;               
r_uranus = 25559;               
r_neptune = 24764;              

// Planet SOI (km)
SOI_mercury = 0.112e6;          
SOI_venus = 0.616e6;            
SOI_earth = 0.924e6;            
SOI_mars = 0.576e6;             
SOI_jupiter = 48.2e6;           
SOI_saturn = 54.6e6;            
SOI_uranus = 51.8e6;            
SOI_neptune = 86.8e6;           



void print_sim_status(struct reb_simulation* const r, float simulation_time){
	printf("----SIMULATION STATUS---- ");

	// print time
	time_t timer;
	char buffer[26];
	struct tm* tm_info;
	time(&timer);
	tm_info = gmtime(&timer);
	strftime(buffer, 27, "%Y-%m-%dT%H:%M:%SZ", tm_info);
	puts(buffer);

	printf("Integration time:               %.8f years\n", simulation_time);
	//printf("Integration time:               %.8f machine time\n", simulation_time*2.0*M_PI);
	printf("Nominal timestep:               %.8f\n", r->dt);
	printf("Current simulation time:        %.8f\n", r->t);
	printf("Number of particles:            %d\n", r->N);
	printf("Number of massive particles:    %d\n", r->N_active);
	printf("Number of test particles:       %d\n", r->N - r->N_active);
	printf("Number of hashed particles:     %d\n", r->hash_ctr);
	printf("Gravity constant:               %.8f\n", r->G);
	//printf("Output data every:              %.8f years\n", r->simulationarchive_interval/2.0/M_PI);
	printf("Output data every:              %.8f years\n", r->simulationarchive_interval);
	printf("Saving to file:                 %s\n",r->simulationarchive_filename);
    
}


// Define our own collision resolve function, which will only record collisions but not change any of the particles.
int collision_record_only(struct reb_simulation* const r, struct reb_collision c){
     double delta_t = 2.*M_PI;
     struct reb_particle* particles = r->particles;
     const double t = r->t;

     // only record a maximum of one collision per year per particle
     if ( particles[c.p1].lastcollision+delta_t < t  &&  particles[c.p2].lastcollision+delta_t < t ){
             particles[c.p1].lastcollision = t;
             particles[c.p2].lastcollision = t;
             printf("\nCollision detected.\n");
             FILE* of = fopen("collisions.txt","a+");                // open file for collision output
             fprintf(of, "%e\t", t);                                 // time
             fprintf(of, "%e\t", (particles[c.p1].x+particles[c.p2].x)/2.);  // x position
             fprintf(of, "%e\t", (particles[c.p1].y+particles[c.p2].y)/2.);  // y position
             fprintf(of, "\n");
             fclose(of);                                             // close file
     }
    return 0;
}


void run_sim_from_archive(char *filename, float simulation_time){
	//char filename[512] = "archive.bin";
	printf("Attempting to read archive file: %s\n",filename);
//	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);
	struct reb_simulation* r = reb_create_simulation_from_binary(filename);
	if (r==NULL){
        	printf("No simulation archive found.\n");
		return;
	}

	printf("Found simulation archive. Loaded snapshot.\n");
	r->simulationarchive_filename = filename;

	/*
	// Setup constants
	r->integrator	= REB_INTEGRATOR_LEAPFROG;
	r->gravity	= REB_GRAVITY_BASIC;
	r->boundary	= REB_BOUNDARY_OPEN;
	r->opening_angle2	= 1.5;	// This constant determines the accuracy of the tree code gravity estimate.
	r->G 		= 1;		
	r->softening 	= 0.02;		// Gravitational softening length
	r->dt 		= 3e-2;		// Timestep
	const double boxsize = 10.2;
	reb_configure_box(r,boxsize,1,1,1);
	*/
	
	print_sim_status(r, simulation_time);
	
	double original_dt = r->dt;

	// set up low step integration to get out of Earth's SOI
	const double earth_SOI = 0.006177; // AU
	const int earth = 3;
	const int meteoroid = r->N_active;
	double prx  = r->particles[meteoroid].x - r->particles[earth].x;
	double pry  = r->particles[meteoroid].y - r->particles[earth].y;
	double prz  = r->particles[meteoroid].z - r->particles[earth].z;
	double dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	printf("Decreasing timestep to get out of the Earth SOI.\n");
	while (dist_earth_particle < 3 * earth_SOI){
		r->dt = original_dt / 20.0 * exp(dist_earth_particle/earth_SOI);
		printf("Current dt: %.8f hours - %.8f Earth radii from Earth\n", r->dt*365.25*24.0, dist_earth_particle*1.496e8/6371.0);
    		reb_integrate(r, 50*r->dt);

		double prx  = r->particles[meteoroid].x - r->particles[earth].x;
		double pry  = r->particles[meteoroid].y - r->particles[earth].y;
		double prz  = r->particles[meteoroid].z - r->particles[earth].z;
		dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	}

	printf("Preliminary final encounter resolved.\n");
	print_sim_status(r, simulation_time);

	printf("Switching to main integration\n");

	r->collision            = REB_COLLISION_TREE;
	r->collision_resolve    = collision_record_only;

	r->dt = original_dt;
	reb_integrate(r, simulation_time);
    
	print_sim_status(r, simulation_time);
    
	reb_free_simulation(r);
}

int main(int argc, char* argv[]){
	// Get the number of processors
	// int np = omp_get_num_procs();

	if (argc < 2){
		printf("\n Usage: ./rebound.x ARCHIVEFILE.bin SIMULATION_TIME_EARTH_YEARS\n");
	}

	//int np = atoi(argv[1]);
	char* np_string = getenv("OMP_NUM_THREADS");
        int np = atoi(np_string);
    
	// Set the number of OpenMP threads to be the number of processors	
	omp_set_num_threads(np);

	printf("\nomp_get_num_procs() got %d from env var OMP_NUM_THREADS\n", np);
	
	run_sim_from_archive(argv[1], atof(argv[2]));


}

