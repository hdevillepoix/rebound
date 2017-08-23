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

/*
static const char COLLISION_FILE[]   = "_collisions_record.csv";
static const char ORBIT_CONV_FILE[]  = "_orbits.csv";
static const char ORBIT_ALL_FILE[]   = "_orbits.csv";
*/


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


/**
General function to generate output filename using input file basename
*/
int generate_output_filename(char *filename, struct reb_simulation* const r, const char *output_type){
	strcat(filename, rem_ext(r->simulationarchive_filename));
	
	// cast to int the sim time and use that for the name
	int int_curr_sim_time = (int) r->t;
	char str_curr_sim_time[12];
	sprintf(str_curr_sim_time, "%d", int_curr_sim_time);
	
	if (strcmp(output_type, "COLLISION") == 0){
		strcat(filename, "_collisions_record.csv");
	}
	else if (strcmp(output_type, "OOB") == 0){
		strcat(filename, "_out_of_bounds_orbits.csv");
	}
	else if (strcmp(output_type, "ORBIT_CONVERGENCE") == 0){
		strcat(filename, "_converging_orbits.csv");
	}
	else if (strcmp(output_type, "ORBIT_ALL") == 0){
		//strcat(filename, "_");		
		//strcat(filename, str_curr_sim_time);
		strcat(filename, "_all_orbits.csv");
	}
	else if (strcmp(output_type, "BINARY") == 0){
		strcat(filename, "_");		
		strcat(filename, str_curr_sim_time);
		strcat(filename, ".bin");
	}
	return 0;
}

/**
General function to generate CSV header for output file.
Overwrites existing files
*/
void initialise_output_file(struct reb_simulation* const r, const char *output_type){
	char filename[255] = "";
	generate_output_filename(filename, r, output_type);
	printf("Initialising file: %s\n", filename);

	FILE* of = fopen(filename,"w");
	if (strcmp(output_type, "COLLISION") == 0){
		fprintf(of, "time,massive_particle_hash,test_particle_hash,distance\n");
	}
	else if ((strcmp(output_type, "ORBIT_CONVERGENCE") == 0) || (strcmp(output_type, "ORBIT_ALL") == 0) || (strcmp(output_type, "OOB") == 0)){
		fprintf(of, "time,particle_hash,a,e,i,w,W\n");
		//time, semi-major axis, eccentricity, inclination, Omega (longitude ascending node), omega (argument of pericenter), lambda (mean longitude), period, f (true anomaly).
	}
	fclose(of);
}





/*

int generate_collision_file(char *collision_file, struct reb_simulation* const r){
	//char collision_file[255] = "";
	strcat(collision_file, rem_ext(r->simulationarchive_filename));
	strcat(collision_file, "_collisions_record.csv");
	return 0;
}


int generate_orbit_file(char *orbit_file, struct reb_simulation* const r){
	//char collision_file[255] = "";
	strcat(collision_file, rem_ext(r->simulationarchive_filename));
	strcat(collision_file, "_orbits.csv");
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
*/

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
	particles[0].r = r_sun/AU_KM;

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
int collision_record_and_destroy(struct reb_simulation* const r, struct reb_collision c){
	double delta_t = 1.0; // 1 year
	struct reb_particle* particles = r->particles;
	const double t = r->t;

	// only record if it is not a collision between 2 massive particles (eg. Moon in Earth SOI)
	if (c.p1 >= r->N_active || c.p2 >= r->N_active){

		particles[c.p1].lastcollision = t;
		particles[c.p2].lastcollision = t;
		printf("Collision detected at %f\n", t);
		char collision_file[255] = "";
		generate_output_filename(collision_file, r, "COLLISION");
		//printf("Collision output: %s\n", collision_file);
		FILE* of = fopen(collision_file,"a+");		   // open file for collision output

		int min_index = MIN(c.p1, c.p2);
		int max_index = MAX(c.p1, c.p2);
		// have the massive particle as p1 and the test particle as p2 for easing post mortem analysis
		struct reb_particle diff_p = reb_particle_minus(particles[c.p1],particles[c.p2]);
		double dist_au = sqrt(diff_p.x*diff_p.x+diff_p.y*diff_p.y+diff_p.z*diff_p.z);
		fprintf(of, "%e,%u,%u,%e\n", t, particles[min_index].hash, particles[max_index].hash, dist_au);
		fclose(of);
		// remove particle
		int rmstatus = reb_remove(r, max_index, 1);
		if (rmstatus == 1){
			printf("Particle %d has been removed from simulation after collision with massive particle %d\n", max_index, min_index);
		}
		else{
			printf("PROBLEM REMOVING particle %d after collision with massive particle %d", max_index, min_index);
		}
		
	}
	return 0;
}



/**
Remove particle that are sent on irrealistic orbits
*/
void out_of_bounds_check(struct reb_simulation* const r, struct reb_orbit* orbits){

	char orb_file[255] = "";
	generate_output_filename(orb_file, r, "OOB");
	FILE* of = fopen(orb_file,"a+");

	// iterate all test particles
	// from the end to avoid index problems
	for (int i=r->N-1; i >= r->N_active; i--){
		struct reb_orbit orb = orbits[i];
		if ((orb.e >  0.95) || (orb.a > 6.0) || (orb.a < 0.0)){	
			int hash = r->particles[i].hash;
			fprintf(of, "%e,%u,%e,%e,%e,%e,%e\n", r->t, hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);	
			int rmstatus = reb_remove(r, i, 0);
			if (rmstatus == 1){
				printf("Particle %u has been removed (OOB).\n", hash);
			}
			else{
				printf("PROBLEM REMOVING particle %d\n", i);
			}
		}
	}	
	fclose(of);
}

void source_region_convergence_check(struct reb_simulation* const r, struct reb_orbit* orbits){

	char orb_file[255] = "";
	generate_output_filename(orb_file, r, "ORBIT_CONVERGENCE");
	FILE* of = fopen(orb_file,"a+");

	// iterate all test particles
	for (int i=r->N_active; i < r->N; i++){
		struct reb_orbit orb = orbits[i];
		
		// escape Mars influence: perihelion > Mars aphelion
		if (orb.a*(1-orb.e) > 1.67){
			fprintf(of, "%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);
		}
	}	
	fclose(of);
}


/**w3
DEPRECATED
void source_region_convergence_check(struct reb_simulation* const r){

	char orb_file[255] = "";
	generate_output_filename(orb_file, r, "ORBIT_CONVERGENCE");
	FILE* of = fopen(orb_file,"a+");

	// iterate all test particles
	for (int i=r->N_active; i < r->N; i++){
		struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
		
		// escape Mars influence: perihelion > Mars aphelion
		if (orb.a*(1-orb.e) > 1.67){
			fprintf(of, "%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);
		}
	}	
	fclose(of);
}
*/


/**
DEPRECATED
void save_all_orbits(struct reb_simulation* const r){

	char orb_file[255] = "";
	generate_output_filename(orb_file, r, "ORBIT_ALL");
	FILE* of = fopen(orb_file,"a+");

	// iterate all test particles and massive particles except Sun
	for (int i=1; i < r->N; i++){
		struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
		
		fprintf(of, "%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);		
	}
	fclose(of);
}
*/
void save_all_orbits(struct reb_simulation* const r, struct reb_orbit* orbits){

	char orb_file[255] = "";
	generate_output_filename(orb_file, r, "ORBIT_ALL");
	FILE* of = fopen(orb_file,"a+");

	// iterate all test particles and massive particles except Sun
	for (int i=1; i < r->N; i++){
		struct reb_orbit orb = orbits[i];
		
		fprintf(of, "%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);		
	}
	fclose(of);
}


/**
All functions that require orbital elements to be calculated
*/
void orbital_check(struct reb_simulation* const r){
	
	// calculate orbits of all particles
	struct reb_orbit* orbits = malloc(r->N * sizeof(struct reb_orbit));
	
	#pragma omp parallel for
	for (int i=1; i < r->N; i++){
		orbits[i] = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);	
	}
	
	// record orbits of all surviving particles
	save_all_orbits(r, orbits);

	// check if some particles have attained convergence to MB
	source_region_convergence_check(r, orbits);
	
	
	// remove particles that have gone outside of the plausible orbital space, after initial close encounter
	if (r->t < -0.5){
		out_of_bounds_check(r, orbits);
	}
	
	free(orbits);
}
		





int simulation_health_check(struct reb_simulation* const r){
	// check that there are still test particles in the simulation
	// exit if not
	
	if (r->N_active == r->N){
		printf("Simulation has ran out of test particles.\nExiting at t = %.8f years ...\n", r->t);
		exit(2);
	}
	return 0;
	
	
	
}

void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, 20.)){
		//reb_output_timing(r, 100);
		
		
		// orbital check stuff
		orbital_check(r);
		
		// only snapshot every 1000 years
		if (((int) r->t) % 1000 == 0){
			char binary_file[255] = "";
			generate_output_filename(binary_file, r, "BINARY");
			reb_output_binary(r, binary_file);
			printf("t= %f years. N_active/N = %d/%d. Snapshot: %s\n", r->t, r->N_active, r->N, binary_file);
		}
		
		// kill simulation if all test particles have died
		simulation_health_check(r);
		//reb_output_orbits(r, "orbits.csv");
		// save a snapshot of the simulation in case it gets interrupted
		//reb_output_binary(struct reb_simulation *r, char *filename)
	}
}


void run_sim_from_archive_start(char *filename, float simulation_time){
	printf("Initiating pre-set simulation.\nReading archive file: %S\n",filename);
	struct reb_simulation* r = reb_create_simulation_from_simulationarchive(filename);
//	struct reb_simulation* r = reb_create_simulation_from_binary(filename);
	if (r==NULL){
        	printf("No simulation archive found.\n");
		exit(1);
	}

	printf("Found simulation archive. Loaded pre-made simulation set up.\n");
	
	// set to give output functions access to basename
	r->simulationarchive_filename = filename;

	initialise_output_file(r, "COLLISION");
	initialise_output_file(r, "ORBIT_CONVERGENCE");
	initialise_output_file(r, "ORBIT_ALL");
	initialise_output_file(r, "OOB");
	


	//r->collision            = REB_COLLISION_TREE;
	//r->integrator           = REB_INTEGRATOR_IAS15;
	//r->integrator           = REB_INTEGRATOR_WHFAST;
	
	
	// assignes uint hashes to all particles
	for (int i=0; i < r->N; i++){
		r->particles[i].hash = i;
	}
	
	//printf("Integrator: %s\n",r->integrator);
	//r->integrator.min_dt = 1./24./365.25; //1 hour
	r->heartbeat            = heartbeat;
	r->simulationarchive_interval = -1.0;
	r->simulationarchive_interval = 10e10;


	print_sim_status(r, simulation_time);
	

	//const double original_dt = r->dt;
	const double original_dt = -0.01; // TODO FIXME
	
	
	printf("Switching to IAS15 integrator for Earth approach resolution.\n");
	r->integrator           = REB_INTEGRATOR_IAS15;

	// set up low step integration to get out of Earth's SOI
	const double SOI_earth_AU = SOI_earth/AU_KM; // AU
	const int earth = 3;
	const int meteoroid = r->N_active+1;
	double prx  = r->particles[meteoroid].x - r->particles[earth].x;
	double pry  = r->particles[meteoroid].y - r->particles[earth].y;
	double prz  = r->particles[meteoroid].z - r->particles[earth].z;
	double dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	
	while (dist_earth_particle < 3 * SOI_earth_AU){
		r->dt = original_dt / 20.0 * exp(dist_earth_particle/SOI_earth_AU);
		printf("Current dt: %.8f hours - %.8f Earth radii from Earth\n", r->dt*365.25*24.0, dist_earth_particle*1.496e8/6371.0);
    		reb_integrate(r, r->t - 0.001);

		double prx  = r->particles[meteoroid].x - r->particles[earth].x;
		double pry  = r->particles[meteoroid].y - r->particles[earth].y;
		double prz  = r->particles[meteoroid].z - r->particles[earth].z;
		dist_earth_particle = sqrt(prx*prx + pry*pry + prz*prz);
	}
	
	// save to binary after the resolution of the Earth encounter
	char binary_file[255] = "";
	generate_output_filename(binary_file, r, "BINARY");
	reb_output_binary(r, binary_file);
	
	// save all orbits after the resolution of the Earth encounter
		
		
	
	printf("Preliminary final encounter resolved.\n");
	print_sim_status(r, simulation_time);
	
	
	printf("Switching to WHFAST integrator for main integration.\n");
	r->integrator           = REB_INTEGRATOR_WHFAST;
	r->dt = original_dt;
	
	/*
	printf("Setting hill sphere collisions.\n");
	//r->collision            = REB_COLLISION_DIRECT;
	r->collision            = REB_COLLISION_TREE;
	r->collision_resolve    = collision_record_and_destroy;
	increase_radius_for_close_encounter_record(r);
	
	// check virtual radii of massive particles
	for (int i=0; i < r->N_active; i++){
		printf("Particle %d has a virtual radius of %.8f km\n", i, r->particles[i].r*AU_KM);
	}
	
	*/
	r->collision            = REB_COLLISION_NONE;

	reb_integrate(r, -simulation_time);
    
	print_sim_status(r, simulation_time);
    
	reb_free_simulation(r);
}



void run_sim_from_archive_snapshot(char *filename, float simulation_time, char *snapshot_fname){
	printf("Initiating simulation restart.\nReading archive file: %s\n",snapshot_fname);
	struct reb_simulation* r = reb_create_simulation_from_binary(snapshot_fname);
	if (r==NULL){
        	printf("No simulation archive found.\n");
		exit(1);
	}

	printf("Found simulation archive. Loaded snapshot at %f years.\n", r->t);
	
	/*
	for (int i=1; i < r->N; i++){
		struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
		printf("%e,%u,%e,%e,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, r->particles[i].x, r->particles[i].y, r->particles[i].z, r->particles[i].vx, r->particles[i].vy, r->particles[i].vz);
		printf("%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);
	}	
	*/
	
	// keep as basename
	r->simulationarchive_filename = filename;

	r->heartbeat            = heartbeat;
	r->simulationarchive_interval = -1.0;
	r->simulationarchive_interval = 10e10;


	print_sim_status(r, simulation_time);
	
	
	printf("Switching to WHFAST integrator for main integration.\n");
	r->integrator           = REB_INTEGRATOR_WHFAST;
	
	/*
	printf("Setting hill sphere collisions.\n");
	//r->collision            = REB_COLLISION_DIRECT;
	r->collision            = REB_COLLISION_TREE;
	r->collision_resolve    = collision_record_and_destroy;
	increase_radius_for_close_encounter_record(r);
	
	// check virtual radii of massive particles
	for (int i=0; i < r->N_active; i++){
		printf("Particle %d has a virtual radius of %.8f km\n", i, r->particles[i].r*AU_KM);
	}
	
	*/
	r->collision            = REB_COLLISION_NONE;



	reb_integrate(r, -simulation_time);
    
	print_sim_status(r, simulation_time);
	
	for (int i=r->N_active; i < r->N; i++){
		struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
		printf("%e,%u,%e,%e,%e,%e,%e\n", r->t, r->particles[i].hash, orb.a, orb.e, orb.inc, orb.omega, orb.Omega);
	}
    
	reb_free_simulation(r);
}


int main(int argc, char* argv[]){
	// Get the number of processors
	// int np = omp_get_num_procs();

	if (argc < 3){
		printf("\n Usage: ./rebound.x ARCHIVEFILE.bin SIMULATION_TIME_EARTH_YEARS [snapshot.bin]\n");
		exit(1);
	}

	//int np = atoi(argv[1]);
	char* np_string = getenv("OMP_NUM_THREADS");
        int np = atoi(np_string);
    
	// Set the number of OpenMP threads to be the number of processors	
	omp_set_num_threads(np);

	printf("\nomp_get_num_procs() got %d from env var OMP_NUM_THREADS\n", np);
	
	printf("Got %d arguments:\n", argc);
	for (int i = 0; i < argc; i++) {
		printf("- argument %d: %s\n", i, argv[i]);
	}
	
	if (argc == 3){
		printf("\n Usage: ./rebound.x ARCHIVEFILE.bin SIMULATION_TIME_EARTH_YEARS\n");
		run_sim_from_archive_start(argv[1], atof(argv[2]));
	}
	else if (argc == 4){
		run_sim_from_archive_snapshot(argv[1], atof(argv[2]), argv[3]);
	}
	else{
		printf("\nInvalid number of arguments\n");
		exit(1);
	}
	


}

