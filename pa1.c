
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <getopt.h>
#include <omp.h>
#include <sqlite3.h>

#define FILEHEADER "index, mass, x, y, z, vx, vy, vz\n"

struct particle {
	double mass;
	double pos[2];
	double vel[2];
};

/* Initializes a set of particles which are uniformly distributed over a unit cube */
struct particle *initParticles(int pnum);
struct particle *loadParticles(const char *dbname, int *pnum);
void saveParticles(struct particle *parts, int pnum, const char *fname);
void slowersim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
void slowsim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
double calcMomentum(struct particle *particles, int pnum);
double calcEnergy(struct particle *particles, int pnum, double gravity);

double magsq(double x, double y);

int runSim(int pnum, int maxstep, double timestep, double gravity,
					 bool save, bool load)
{
	struct particle *particles = NULL;
	if(!load) {
		particles = initParticles(pnum);
		if(save)
			saveParticles(particles, pnum, "start.db");
	}
	else {
		particles = loadParticles("start.db", &pnum);
	}
	struct particle *copy = (struct particle *)malloc(sizeof(struct particle[pnum]));
	memcpy(copy, particles, sizeof(struct particle[pnum]));
	printf("%5d, %4d, %10.9f, %10.9f, %10.9lf, %10.9lf, ",
				 pnum, maxstep, timestep, gravity,
				 calcMomentum(particles, pnum),
				 calcEnergy(particles, pnum, gravity));
	struct timeval tv1;
	gettimeofday(&tv1, NULL);
	slowersim(particles, pnum, maxstep, timestep, gravity);
	struct timeval tv2, elapsed;
	gettimeofday(&tv2, NULL);
	timersub(&tv2, &tv1, &elapsed);
	printf(" %u.%06u, ", elapsed.tv_sec, elapsed.tv_usec);
	printf("%10.9lf, %10.9lf, ",
				 calcMomentum(particles, pnum),
				 calcEnergy(particles, pnum, gravity));
	free(particles);

	gettimeofday(&tv1, NULL);
	particles = copy;
	slowsim(particles, pnum, maxstep, timestep, gravity);
	gettimeofday(&tv2, NULL);
	timersub(&tv2, &tv1, &elapsed);
	printf("%u.%06u, ", elapsed.tv_sec, elapsed.tv_usec);
	printf("%10.9lf, %10.9lf\n",
				 calcMomentum(particles, pnum),
				 calcEnergy(particles, pnum, gravity));
	free(particles);
	return 0;
}

void slowersim(struct particle *particles, int pnum, int maxstep,
							 double timestep, double gravity)
{
	register double ax, ay;
	for(int s = 0; s < maxstep; s++) {
		for(int i = 0; i < pnum; i++) {
			ax = ay = 0.0;
			register double x = particles[i].pos[0],
				y = particles[i].pos[1];
			for(int j = 0; j < pnum; j++) {
				if(i == j)
					continue;
				register double dispx = x - particles[j].pos[0],
					dispy = y - particles[j].pos[1];
				double rsq = dispx * dispx + dispy * dispy,
					c = gravity * particles[j].mass / (rsq * sqrt(rsq));
				ax -= c * dispx;
				ay -= c * dispy;
			}
			particles[i].vel[0] += ax * timestep;
			particles[i].vel[1] += ay * timestep;
		}
		for(int i = 0; i < pnum; i++) {
			particles[i].pos[0] += particles[i].vel[0] * timestep;
			particles[i].pos[1] += particles[i].vel[1] * timestep;
		}
	}
}

void slowsim(struct particle *particles, int pnum, int maxstep,
							 double timestep, double gravity)
{
	double *force = malloc(sizeof(double[pnum * 2]));
	for(int s = 0; s < maxstep; s++) {
		for(int i = 0; i < pnum; i++) {
			for(int k = 0; k < 2; k++)
				force[2 * i + k] = 0;
		}
		for(int i = 0; i < pnum; i++) {
			for(int j = i + 1; j < pnum; j++) {
				register double dispx = particles[i].pos[0] - particles[j].pos[0],
					dispy = particles[i].pos[1] - particles[j].pos[1];
				double c = gravity / ((dispx * dispx + dispy * dispy) *
															sqrt(dispx * dispx + dispy * dispy)),
					a1 = c * particles[j].mass,
					a2 = c * particles[i].mass;
				force[2 * i] -= a1 * dispx;
				force[2 * j + 1] += a2 * dispy;
			}
			particles[i].pos[0] += particles[i].vel[0] * timestep;
			particles[i].pos[1] += particles[i].vel[1] * timestep;
			particles[i].vel[0] += force[2 * i] * timestep;
			particles[i].vel[1] += force[2 * i + 1] * timestep;
		}
	}
	free(force);
}

double calcMomentum(struct particle *parts, int pnum)
{
	double momentum = 0.0;
	for(int i = 0; i < pnum; i++)
		momentum += sqrt(magsq(parts[i].vel[0], parts[i].vel[1])) *
			parts[i].mass;
	return momentum;
}

double calcEnergy(struct particle *parts, int pnum, double gravity)
{
	double energy = 0;
	for(int i = 0; i < pnum; i++) {
		energy += magsq(parts[i].vel[0], parts[i].vel[1]) *
			parts[i].mass / 2;
		for(int j = i + 1; j < pnum; j++) {
			double tmpmass = parts[i].mass * parts[j].mass,
				tmpdist = sqrt(magsq(parts[i].pos[0] - parts[j].pos[0],
														 parts[i].pos[1] - parts[j].pos[1]));
			energy += -gravity * tmpmass / tmpdist;
		}
	}
	return energy;
}

double magsq(double x, double y)
{
	return x * x + y * y;
}

struct particle *initParticles(int pnum)
{
	struct particle *particles = (struct particle *)
		malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
    particles[i].mass = ((double)rand()) / RAND_MAX;
		particles[i].pos[0] = ((double)rand() - 1) / RAND_MAX;
		particles[i].pos[1] = ((double)rand() - 1) / RAND_MAX;
		particles[i].vel[0] = ((double)rand() - 1) / RAND_MAX;
		particles[i].vel[1] = ((double)rand() - 1) / RAND_MAX;
	}
	return particles;
}

struct particle *loadParticles(const char *fname, int *pnum)
{
	sqlite3 *dbh;
	int result = sqlite3_open(fname, &dbh);
	if(result != SQLITE_OK) {
		printf("Error: loadParticles: Could not open database!\n");
		sqlite3_close(dbh);
		return NULL;
	}
	sqlite3_stmt *stmt;
	const char getCount[] = "SELECT COUNT(*) FROM Particles;";
	result = sqlite3_prepare_v2(dbh, getCount, sizeof(getCount), &stmt, NULL);
	if(result != SQLITE_OK) {
		printf("Error: loadParticles: Could count the number of rows!\n");
		sqlite3_close(dbh);
		return NULL;
	}
	sqlite3_step(stmt);
	*pnum = sqlite3_column_int(stmt, 0);
	sqlite3_finalize(stmt);
	if(*pnum == 0) {
		sqlite3_close(dbh);
		return NULL;
	}

	struct particle *particles = malloc(sizeof(struct particle[*pnum]));
	const char getRows[] = "SELECT * FROM Particles;";
	result = sqlite3_prepare_v2(dbh, getRows, sizeof(getRows), &stmt, NULL);
	if(result != SQLITE_OK) {
		printf("Error: loadParticles: Could get the rows!\n");
		sqlite3_close(dbh);
		return NULL;
	}
	result = sqlite3_step(stmt);
	for(int i = 0; i < *pnum && result == SQLITE_ROW; i++) {
		particles[i].mass = sqlite3_column_double(stmt, 1);
		particles[i].pos[0] = sqlite3_column_double(stmt, 2);
		particles[i].pos[1] = sqlite3_column_double(stmt, 3);
		particles[i].vel[0] = sqlite3_column_double(stmt, 4);
		particles[i].vel[1] = sqlite3_column_double(stmt, 5);
		
		result = sqlite3_step(stmt);
	}
	sqlite3_finalize(stmt);
	sqlite3_close(dbh);
	return particles;
}

void saveParticles(struct particle *parts, int pnum, const char *fname)
{
	sqlite3 *dbh;
	int err = sqlite3_open(fname, &dbh);
	if(err != SQLITE_OK) {
		printf("Error: saveParticles: Could not open database!\n");
		sqlite3_close(dbh);
		return;
	}
	const char dropTable[] = "DROP TABLE IF EXISTS Particles;";
	sqlite3_stmt *stmt;
	err = sqlite3_prepare_v2(dbh, dropTable, sizeof(dropTable), &stmt, NULL);
	if(err != SQLITE_OK) {
		printf("Error: saveParticles: Could not create statement!\n");
		sqlite3_close(dbh);
		return;
	}
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	const char createTable[] = "CREATE TABLE IF NOT EXISTS Particles"
		"(id Integer, mass Real, x Real, y Real, vx Real, vy Real);";
	err = sqlite3_prepare_v2(dbh, createTable, sizeof(createTable), &stmt, NULL);
	if(err != SQLITE_OK) {
		printf("Error: saveParticles: Could not create statement!\n");
		sqlite3_close(dbh);
		return;
	}
	sqlite3_step(stmt);
	sqlite3_finalize(stmt);
	const char *insert = "INSERT INTO Particles VALUES (%d, %.20f, %.20f, %.20f, %.20f, %.20f);";
	const size_t buflen = sizeof(char[256]);
	char *query = (char *)malloc(buflen);
	for(int i = 0; i < pnum; i++) {
		snprintf(query, buflen, insert, i, parts[i].mass,
						 parts[i].pos[0], parts[i].pos[1],
						 parts[i].vel[0], parts[i].vel[1]);
		query[buflen - 1] = 0;
		err = sqlite3_prepare_v2(dbh, query, buflen, &stmt, NULL);
		if(err != SQLITE_OK) {
			printf("Error: saveParticles: Could not create insert statement!\n");
			sqlite3_close(dbh);
			return;
		}
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	}
	free(query);
	sqlite3_close(dbh);
}

int main(int argc, char **argv)
{
	/* Default Values */
	int pnum = 1000, maxstep = 4;
	double ts = 0.0001, gravity = 0.0000000000667384;
	bool load = false, save = false;
	for(;;) {
		static struct option longOpts[] =
			{
				{"count", required_argument, 0, 'c'},
				{"timestep", required_argument, 0, 't'},
				{"steps", required_argument, 0, 's'},
				{"gravity", required_argument, 0, 'g'},
				{"load", no_argument, 0, 'l'},
				{"save", no_argument, 0, 'w'},
				{0, 0, 0, 0}
			};
		int index = 0;
		int curopt = getopt_long(argc, argv, "c:t:s:g:", longOpts, &index);
		if(curopt == -1)
			break;
		switch(curopt) {
		case 'c':
			sscanf(optarg, "%d", &pnum);
			break;
		case 't':
			sscanf(optarg, "%lf", &ts);
			break;
		case 's':
			sscanf(optarg, "%d", &maxstep);
			break;
		case 'g':
			sscanf(optarg, "%lf", &gravity);
			break;
		case 'l':
			load = true;
			break;
		case 'w':
			save = true;
		}
	}
	if(pnum == 0) {
		printf("Cannot simulate less than 1 particle. Exiting...");
		return 1;
	}
	if(maxstep == 0) {
		printf("Cannot simulate less than 1 step. Exiting...");
		return 1;
	}
	if(ts == 0.0) {
		printf("Cannot simulate with a zero timestep. Exiting...");
		return 1;
	}
	srand(time(0));
	return runSim(pnum, maxstep, ts, gravity, save, load);
}
