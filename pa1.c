
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
	double pos[2];
	double vel[2];
};

struct particle *initParticles(int pnum);
void saveParticles(struct particle *parts, int pnum, const char *fname);
void slowersim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
void slowsim(struct particle *particles, int pnum, int maxstep,
						 double timestep, double gravity);
double calcMomentum(struct particle *particles, int pnum);
double calcEnergy(struct particle *particles, int pnum, double gravity);

double magsq(double x, double y)
{
	return x * x + y * y;
}

int runSim(int pnum, int maxstep, double timestep, double gravity)
{
	struct particle *particles = initParticles(pnum),
		*copy = (struct particle *)malloc(sizeof(struct particle[pnum]));
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

	particles = copy;
	slowsim(particles, pnum, maxstep, timestep, gravity);
	gettimeofday(&tv1, NULL);
	timersub(&tv1, &tv2, &elapsed);
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
	register double fx, fy;
	for(int s = 0; s < maxstep; s++) {
		for(int i = 0; i < pnum; i++) {
			fx = fy = 0.0;
			register double x = particles[i].pos[0],
				y = particles[i].pos[1];
			for(int j = 0; j < pnum; j++) {
				if(i == j)
					continue;
				double r = sqrt(x * x + y * y);
				register double c = gravity * particles[j].mass / (r * r * r);
				fx -= c * (x - particles[j].pos[0]);
				fy -= c * (y - particles[j].pos[1]);
			}
			particles[i].vel[0] += fx * timestep;
			particles[i].vel[1] += fy * timestep;
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
				double r = sqrt(magsq(particles[i].pos[0] - particles[j].pos[0],
															particles[i].pos[1] - particles[j].pos[1]));
				double c = gravity / (r * r * r),
					a1 = c * particles[j].mass,
					a2 = c * particles[i].mass;
				for(int k = 0; k < 2; k++) {
					force[2 * i + k] -= a1 * (particles[i].pos[k] - particles[j].pos[k]);
					force[2 * j + k] += a2 * (particles[i].pos[k] - particles[j].pos[k]);
				}
			}
			for(int k = 0; k < 2; k++) {
				double dvdt = force[2 * i + k];
				particles[i].pos[k] += (particles[i].vel[k] + dvdt / 2) * timestep;
				particles[i].vel[k] += dvdt * timestep;
			}
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

struct particle *loadParticles(char *fname)
{
	FILE *file = fopen(fname, "r");
	int pnum = 0;
	fscanf(file, "%d", &pnum);
	fseek(file, SEEK_CUR, sizeof(FILEHEADER));
	struct particle *particles = (struct particle *)
		malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
		int tmp;
		fscanf(file, "%4d, %10lf, "
					 "%10lf, %10lf, "
					 "%10lf, %10lf\n",
					 &tmp, &particles[i].mass,
					 &particles[i].pos[0], &particles[i].pos[1],
					 &particles[i].vel[0], &particles[i].vel[1]);
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
						"%10.9f, %10.9f, "
						"%10.9f, %10.9f\n",
						i, parts[i].mass,
						parts[i].pos[0], parts[i].pos[1],
						parts[i].vel[0], parts[i].vel[1]);
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
