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


void run_sim_from_archive(char *filename, float simulation_time){
	//char filename[512] = "archive.bin";

//	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);
	struct reb_simulation* r = reb_create_simulation_from_binary(filename);
	if (r==NULL){
        	printf("No simulation archive found.\n");
		return;
	}

	printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
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
	

	printf("Integration time: %.8f years\n",simulation_time);
	printf("Nominal timestep: %.8f\n",r->dt);
	printf("Number of particles: %d\n",r->N);
        printf("Number of massive particles: %d\n",r->N_active);
	printf("Output data every %.8f years\n",r->simulationarchive_interval/2.0/M_PI);
        printf("Integration time: %.8f (machine time)\n",simulation_time*2.0*M_PI);
	printf("Saving to file: %s\n",r->simulationarchive_filename);

	reb_integrate(r, simulation_time*2.0*M_PI);
	reb_free_simulation(r);
}

int main(int argc, char* argv[]){
	// Get the number of processors
	// int np = omp_get_num_procs();

	if (argc < 3){
		printf("\n Usage: ./rebound.x OMP_NUM_PROCS ARCHIVEFILE.bin SIMULATION_TIME\n");
	}

	int np = atoi(argv[1]);
	// Set the number of OpenMP threads to be the number of processors	
	omp_set_num_threads(np);

	printf("\n omp_get_num_procs() = %d \n", np);
	
	run_sim_from_archive(argv[2], atof(argv[3]));


}

