/*
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
#include <string.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// 1 AU in km
const float AU_KM = 149597870.7;

// Planet Radii (km)
const float r_sun = 695700.0;  
const float r_mercury = 2439.0;		  
const float r_venus = 6052.0;		    
const float r_earth = 6378.1363;     
const float r_moon = 1737.0;       
const float r_mars = 3397.3;		   
const float r_jupiter = 71492.0;		 
const float r_saturn = 60268.0;		  
const float r_uranus = 25559.0;		  
const float r_neptune = 24764.0;		 
 
// Planet SOI (km)
const float SOI_mercury = 0.112e6;          
const float SOI_venus = 0.616e6;            
const float SOI_earth = 0.924e6;  
const float SOI_moon = 0.0661e6;
const float SOI_mars = 0.576e6;		
const float SOI_jupiter = 48.2e6;           
const float SOI_saturn = 54.6e6;            
const float SOI_uranus = 51.8e6;            
const float SOI_neptune = 86.8e6;    




void print_sim_status(struct reb_simulation* const r, float simulation_time){
	printf("\n---- SIMULATION STATUS ---- ");

	// print time
	time_t timer;
	char buffer[26];
	struct tm* tm_info;
	time(&timer);
	tm_info = gmtime(&timer);
	strftime(buffer, 27, "%Y-%m-%dT%H:%M:%SZ", tm_info);
	puts(buffer);

	printf("Integration time:		%.8f years\n", simulation_time);
	//printf("Integration time:		  %.8f machine time\n", simulation_time*2.0*M_PI);
	printf("Nominal timestep:		%.8f\n", r->dt);
	printf("Current simulation time:        %.8f\n", r->t);
	printf("Number of particles:            %d\n", r->N);
	printf("Number of massive particles:    %d\n", r->N_active);
	printf("Number of test particles:       %d\n", r->N - r->N_active);
	printf("Number of hashed particles:     %d\n", r->hash_ctr);
	printf("Gravity constant:		%.8f\n", r->G);
	//printf("Output data every:		  %.8f years\n", r->simulationarchive_interval/2.0/M_PI);
	printf("Output data every:		%.8f years\n", r->simulationarchive_interval);
	printf("Saving to file:		        %s\n\n",r->simulationarchive_filename);
    
}

/*
Remove extension from path

*/
char *rem_ext(char* mystr) {
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}


int generate_collision_file(char *collision_file, struct reb_simulation* const r){
	//char collision_file[255] = "";
	strcat(collision_file, rem_ext(r->simulationarchive_filename));
	strcat(collision_file, "_collisions_record.csv");
	return 0;
}

void initialise_collision_file(struct reb_simulation* const r){
	char collision_file[255] = "";
	generate_collision_file(collision_file, r);
	printf("Initialising collision file: %s\n", collision_file);

	FILE* of = fopen(collision_file,"w");		   // open file for collision output
	fprintf(of, "time,particle_1,particle_2,distance\n");				       // time
	fclose(of);
}


/*
    sim.add("Sun", date=date)
    sim.add("Mercury", date=date)
    sim.add("Venus", date=date)
    sim.add("399", date=date) # Earth
    sim.add("301", date=date) # Moon
    sim.add("Mars", date=date)
    sim.add("Jupiter", date=date)
    sim.add("Saturn", date=date)
    sim.add("Uranus", date=date)
    sim.add("Neptune", date=date)       
*/



int increase_radius_for_close_encounter_record(struct reb_simulation* r){
	
	struct reb_particle* particles = r->particles;
	particles[1].r = SOI_mercury/AU_KM;
	particles[2].r = SOI_venus/AU_KM;
	particles[3].r = SOI_earth/AU_KM;
	particles[4].r = SOI_moon/AU_KM;
	particles[5].r = SOI_mars/AU_KM;
	particles[6].r = SOI_jupiter/AU_KM;
	particles[7].r = SOI_saturn/AU_KM;
	particles[8].r = SOI_uranus/AU_KM;
	particles[9].r = SOI_neptune/AU_KM;

	
	//particles[15].r = 1.0/AU_KM;

	printf("Planets now have their SOI as radius.\n");

	return 0;
}

// Define our own collision resolve function, which will only record collisions but not change any of the particles.
int collision_record_only(struct reb_simulation* const r, struct reb_collision c){
	double delta_t = 1.0; // 1 year
	struct reb_particle* particles = r->particles;
	const double t = r->t;

	// only record a maximum of one collision per year per particle and if it is not a collision between 2 massive particles (eg. Moon in Earth SOI)
	if ( (c.p1 >= r->N_active || c.p2 >= r->N_active) && particles[c.p1].lastcollision+delta_t < t  &&  particles[c.p2].lastcollision+delta_t < t ){
		particles[c.p1].lastcollision = t;
		particles[c.p2].lastcollision = t;
		printf("Collision detected at %f\n", t);
		char collision_file[255] = "";
		generate_collision_file(collision_file, r);
		//printf("Collision output: %s\n", collision_file);
		FILE* of = fopen(collision_file,"a+");		   // open file for collision output

		int min_index = MIN(c.p1, c.p2);
		int max_index = MAX(c.p1, c.p2);
		// have the massive particle as p1 and the test particle as p2 for easing post mortem analysis
		struct reb_particle diff_p = reb_particle_minus(particles[c.p1],particles[c.p2]);
		fprintf(of, "%e,%d,%d,%e\n", t, min_index, max_index,sqrt(diff_p.x*diff_p.x+diff_p.y*diff_p.y+diff_p.z*diff_p.z));
		fclose(of);
	}
	return 0;
}


void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, 0.1)){
		reb_output_timing(r, 100);
		//reb_output_orbits(r, "orbits.csv");
	}
}


void run_sim_from_archive(char *filename, float simulation_time){
	//char filename[512] = "archive.bin";
	printf("Attempting to read archive file: %s\n",filename);
	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);
//	struct reb_simulation* r = reb_create_simulation_from_binary(filename);
	if (r==NULL){
        	printf("No simulation archive found.\n");
		return;
	}

	printf("Found simulation archive. Loaded snapshot.\n");
	r->simulationarchive_filename = filename;


	initialise_collision_file(r);

	//r->collision            = REB_COLLISION_TREE;
	//r->integrator           = REB_INTEGRATOR_IAS15;
	//r->integrator           = REB_INTEGRATOR_WHFAST;
	r->collision            = REB_COLLISION_DIRECT;
	r->collision_resolve    = collision_record_only;
	increase_radius_for_close_encounter_record(r);
	r->heartbeat            = heartbeat;
	r->simulationarchive_interval = -1.0;
/*
	struct reb_particle* particles = r->particles;
	for (int i=1; i < r->N; i++){
		particles[i].lastcollision = -9.99;
	}
*/


	print_sim_status(r, simulation_time);
	



	const double original_dt = r->dt;

	// set up low step integration to get out of Earth's SOI
	const double SOI_earth_AU = SOI_earth/AU_KM; // AU
	const int earth = 3;
	const int meteoroid = r->N_active+1;
	double prx  = r->particles[meteoroid].x - r->particles[earth].x;
	double pry  = r->particles[meteoroid].y - r->particles[earth].y;
	double prz  = r->particles[meteoroid].z - r->particles[earth].z;
	double dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	printf("Decreasing timestep to get out of the Earth SOI.\n");
	while (dist_earth_particle < 3 * SOI_earth_AU){
		r->dt = original_dt / 20.0 * exp(dist_earth_particle/SOI_earth_AU);
		printf("Current dt: %.8f hours - %.8f Earth radii from Earth\n", r->dt*365.25*24.0, dist_earth_particle*1.496e8/6371.0);
    		reb_integrate(r, -0.001);

		double prx  = r->particles[meteoroid].x - r->particles[earth].x;
		double pry  = r->particles[meteoroid].y - r->particles[earth].y;
		double prz  = r->particles[meteoroid].z - r->particles[earth].z;
		dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	}

	printf("Preliminary final encounter resolved.\n");
	print_sim_status(r, simulation_time);

	printf("Switching to main integration\n");
	


	r->dt = original_dt;



	reb_integrate(r, -simulation_time);
    
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

