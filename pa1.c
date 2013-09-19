
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

#define FILEHEADER "index, mass, x, y, z, vx, vy, vz\n"

struct particle {
	double mass;
	double pos[3];
	double vel[3];
	double force[3];
};

struct particle *initParticles(int pnum);
void saveParticles(struct particle *parts, int pnum, const char *fname);
void slowersim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
void slowsim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
double calcMomentum(struct particle *particles, int pnum);

double magnitude(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

int runSim(int pnum, int maxstep, double timestep, double gravity)
{
	pnum = 2;
	gravity = 1;
	timestep = 0.001;
	printf("Running with pnum: %d, maxstep: %d, timestep: %f, gravity: %f\n",
				 pnum, maxstep, timestep, gravity);
	struct particle *particles = initParticles(pnum),
		*copy = malloc(sizeof(struct particle[pnum]));
	double pos[2][3] = {{0.458650, 0.755605, 0.0}, {0.679296, 0.678865, 0.0}};
	double vel[2][3] = {{-0.649317, -0.478834, 0.0}, {-0.495775, -0.798279, 0.0}};
	double mass[2] = {0.577852, 0.919489};
	for(int i = 0; i < 2; i++) {
		particles[i].mass = mass[i];
		for(int k = 0; k < 3; k++) {
			particles[i].pos[k] = pos[i][k];
			particles[i].vel[k] = vel[i][k];
		}
	}
	memcpy(copy, particles, sizeof(struct particle[pnum]));

	saveParticles(particles, pnum, "Startup.csv");
	printf("Starting Momentum: %10.9lf\n", calcMomentum(particles, pnum));
	struct timeval tv1;
	gettimeofday(&tv1, NULL);
	slowersim(particles, pnum, maxstep, timestep, gravity);
	struct timeval tv2, elapsed;
	gettimeofday(&tv2, NULL);
	timersub(&tv2, &tv1, &elapsed);
	printf("Finished Momentum: %10.9lf\n", calcMomentum(particles, pnum));
	printf("Elapsed time: %u.%06u\n", elapsed.tv_sec, elapsed.tv_usec);
	saveParticles(particles, pnum, "Finish.csv");
	free(particles);

	particles = copy;
	slowsim(particles, pnum, maxstep, timestep, gravity);
	gettimeofday(&tv1, NULL);
	timersub(&tv1, &tv2, &elapsed);
	printf("Finished Momentum: %10.9lf\n", calcMomentum(particles, pnum));
	printf("Elapsed time: %u.%06u\n", elapsed.tv_sec, elapsed.tv_usec);
	saveParticles(particles, pnum, "Finish2.csv");
	free(particles);
	return 0;
}

void slowersim(struct particle *particles, int pnum, int maxstep,
							 double timestep, double gravity)
{
	double time = 0.0;
	for(int s = 0; s < maxstep; s++) {
		printf("Step %d\n", s);
		for(int i = 0; i < pnum; i++) {
			particles[i].force[0] = 0;
			particles[i].force[1] = 0;
			particles[i].force[2] = 0;
			for(int j = 0; j < pnum; j++) {
				if(i == j)
					continue;
				double r = magnitude(particles[i].pos[0] - particles[j].pos[0],
														 particles[i].pos[1] - particles[j].pos[1],
														 particles[i].pos[2] - particles[j].pos[2]);
				for(int k = 0; k < 3; k++)
					particles[i].force[k] += gravity *
						particles[i].mass * particles[j].mass *
						(particles[i].pos[k] - particles[j].pos[k]) / r / r / r;
			}
			for(int k = 0; k < 3; k++) {
				particles[i].pos[k] += (particles[i].vel[k] + particles[i].force[k] /
																2 / particles[i].mass) * timestep;
				particles[i].vel[k] += particles[i].force[k] * timestep
					/ particles[i].mass;
			}
		}
		time += timestep;
	}
}

void slowsim(struct particle *particles, int pnum, int maxstep,
							 double timestep, double gravity)
{
	double time = 0.0;
	for(int s = 0; s < maxstep; s++) {
		for(int i = 0; i < pnum; i++)
			for(int k = 0; k < 3; k++)
				particles[i].force[k] = 0;
		printf("Step %d\n", s);
		for(int i = 0; i < pnum; i++) {
			for(int j = i + 1; j < pnum; j++) {
				double r = magnitude(particles[i].pos[0] - particles[j].pos[0],
														 particles[i].pos[1] - particles[j].pos[1],
														 particles[i].pos[2] - particles[j].pos[2]);
				for(int k = 0; k < 3; k++) {
					double tmp = gravity *
						particles[i].mass * particles[j].mass *
						(particles[i].pos[k] - particles[j].pos[k]) / r / r / r;
					particles[i].force[k] += tmp;
					particles[j].force[k] -= tmp;
				}
			}
			for(int k = 0; k < 3; k++) {
				particles[i].pos[k] += (particles[i].vel[k] + 
					particles[i].force[k] / 2 / particles[i].mass) * timestep;
				particles[i].vel[k] += particles[i].force[k] * timestep /
					particles[i].mass;
			}
		}
		time += timestep;
	}
}

double calcMomentum(struct particle *parts, int pnum)
{
	double momentum = 0.0;
	for(int i = 0; i < pnum; i++)
		momentum += magnitude(parts[i].vel[0], parts[i].vel[1], parts[i].vel[2]) *
			parts[i].mass;
	return momentum;
}

struct particle *initParticles(int pnum)
{
	struct particle *particles = malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
    particles[i].mass = ((double)rand()) / RAND_MAX;
		particles[i].pos[0] = ((double)rand() - 1) / RAND_MAX;
		particles[i].pos[1] = ((double)rand() - 1) / RAND_MAX;
		particles[i].pos[2] = ((double)rand() - 1) / RAND_MAX;
		particles[i].vel[0] = ((double)rand() - 1) / RAND_MAX;
		particles[i].vel[1] = ((double)rand() - 1) / RAND_MAX;
		particles[i].vel[2] = ((double)rand() - 1) / RAND_MAX;
		particles[i].force[0] = ((double)rand() - 1) / RAND_MAX;
		particles[i].force[1] = ((double)rand() - 1) / RAND_MAX;
		particles[i].force[2] = ((double)rand() - 1) / RAND_MAX;
	}
	return particles;
}

struct particle *loadParticles(char *fname)
{
	FILE *file = fopen(fname, "r");
	int pnum = 0;
	fscanf(file, "%d", &pnum);
	fseek(file, SEEK_CUR, sizeof(FILEHEADER));
	struct particle *particles = malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
		int tmp;
		fscanf(file, "%4d, %10lf, "
					 "%10lf, %10lf, %10lf, "
					 "%10lf, %10lf, %10lf\n",
					 &tmp, &particles[i].mass,
					 &particles[i].pos[0], &particles[i].pos[1], &particles[i].pos[2],
					 &particles[i].vel[0], &particles[i].vel[1], &particles[i].vel[2]);
	}
	return particles;
}

void saveParticles(struct particle *parts, int pnum, const char *fname)
{
	FILE *file = fopen(fname, "w");
	fprintf(file, "%d\n", pnum);
	fprintf(file, FILEHEADER);
	for(int i = 0; i < pnum; i++) {
		fprintf(file, "%4d, %10.9f, "
						"%10.9f, %10.9f, %10.9f, "
						"%10.9f, %10.9f, %10.9f\n",
						i, parts[i].mass,
						parts[i].pos[0], parts[i].pos[1], parts[i].pos[2],
						parts[i].vel[0], parts[i].vel[1], parts[i].vel[2]);
	}
	fclose(file);
}

int main(int argc, char **argv)
{
	/* Default Values */
	int pnum = 1000, maxstep = 4;
	double ts = 0.0001, gravity = 0.0000000000667384;
	for(;;) {
		static struct option longOpts[] =
			{
				{"count", required_argument, 0, 'c'},
				{"timestep", required_argument, 0, 't'},
				{"steps", required_argument, 0, 's'},
				{"gravity", required_argument, 0, 'g'},
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
	return runSim(pnum, maxstep, ts, gravity);
}
