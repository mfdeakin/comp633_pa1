
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <getopt.h>

#include "matrix.h"

struct particle {
	struct matrix *pos;
	struct matrix *velocity;
	mtxfp mass;
};

struct particle *initParticles(int pnum);
void saveParticles(struct particle *parts, int pnum, const char *fname);

int runSim(int pnum, int maxstep, mtxfp timestep, mtxfp gravity)
{
	printf("Running with pnum: %d, maxstep: %d, timestep: %f, gravity: %f\n",
				 pnum, maxstep, timestep, gravity);
	struct particle *particles = initParticles(pnum);
	
	saveParticles(particles, pnum, "Startup.csv");
	mtxfp time = 0;
	struct timeval tv1, tv2, elapsed;
	gettimeofday(&tv1, NULL);
	//Leak not in this loop
	for(int s = 0; s < maxstep; s++) {
		printf("Step %d\n", s);
		for(int i = 0; i < pnum; i++) {
			struct matrix *force = mtxCreate(1, 3);
			for(int j = 0; j < pnum; j++) {
				if(i == j)
					continue;
				struct matrix *disp = mtxSub(particles[i].pos, particles[j].pos);
				mtxfp xdisp = mtxGet(disp, 0, 0),
					ydisp = mtxGet(disp, 0, 1),
					zdisp = mtxGet(disp, 0, 2),
					dist = sqrt(xdisp * xdisp + ydisp * ydisp + zdisp * zdisp);
				struct matrix *fij = mtxScalarMul(disp, -gravity * particles[i].mass *
																					particles[j].mass / dist / dist / dist),
					*newforce = mtxAdd(force, fij);
				mtxFree(disp);
				mtxFree(fij);
				mtxFree(force);
				force = newforce;
			}
			struct matrix *dvdt = mtxScalarMul(force, timestep / particles[i].mass),
				*newvel = mtxAdd(particles[i].velocity, dvdt),
				*dxdt = mtxScalarMul(newvel, timestep),
				*newpos = mtxAdd(particles[i].pos, dxdt);
			mtxFree(force);
			mtxFree(dvdt);
			mtxFree(dxdt);
			mtxFree(particles[i].velocity);
			mtxFree(particles[i].pos);
			particles[i].velocity = newvel;
			particles[i].pos = newpos;
		}
		time += timestep;
	}
	gettimeofday(&tv2, NULL);
	timersub(&tv2, &tv1, &elapsed);
	printf("Elapsed time: %u.%06u\n", elapsed.tv_sec, elapsed.tv_usec);
	saveParticles(particles, pnum, "Finish.csv");
	
	for(int i = 0; i < pnum; i++) {
		mtxFree(particles[i].pos);
		mtxFree(particles[i].velocity);
	}
	free(particles);
	return 0;
}

struct particle *initParticles(int pnum)
{
	struct particle *particles = malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
		particles[i].pos = mtxCreate(1, 3);
		particles[i].velocity = mtxCreate(1, 3);
		mtxSet(particles[i].pos, 0, 0, ((mtxfp)(rand() - 1)) / RAND_MAX);
		mtxSet(particles[i].pos, 0, 1, ((mtxfp)(rand() - 1)) / RAND_MAX);
		mtxSet(particles[i].pos, 0, 2, ((mtxfp)(rand() - 1)) / RAND_MAX);
		mtxSet(particles[i].velocity, 0, 0, ((mtxfp)(rand() - 1)) / RAND_MAX);
		mtxSet(particles[i].velocity, 0, 1, ((mtxfp)(rand() - 1)) / RAND_MAX);
		mtxSet(particles[i].velocity, 0, 2, ((mtxfp)(rand() - 1)) / RAND_MAX);
    particles[i].mass = ((mtxfp)rand()) / RAND_MAX;
	}
	return particles;
}

void saveParticles(struct particle *parts, int pnum, const char *fname)
{
	FILE *file = fopen(fname, "w");
	fprintf(file, "index, mass, x, y, z, vx, vy, vz\n");
	for(int i = 0; i < pnum; i++) {
		fprintf(file, "%4d, %10.9f, "
						"%10.9f, %10.9f, %10.9f, "
						"%10.9f, %10.9f, %10.9f\n",
						i, parts[i].mass,
						mtxGet(parts[i].pos, 0, 0),
						mtxGet(parts[i].pos, 0, 1),
						mtxGet(parts[i].pos, 0, 2),
						mtxGet(parts[i].velocity, 0, 0),
						mtxGet(parts[i].velocity, 0, 1),
						mtxGet(parts[i].velocity, 0, 2));
	}
	fclose(file);
}

int main(int argc, char **argv)
{
	/* Default Values */
	int pnum = 1000, maxstep = 4;
	mtxfp ts = 0.0001, gravity = 0.0000000000667384;
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
		double tmp;
		switch(curopt) {
		case 'c':
			sscanf(optarg, "%d", &pnum);
			break;
		case 't':
			sscanf(optarg, "%lf", &tmp);
			ts = (mtxfp)tmp;
			break;
		case 's':
			sscanf(optarg, "%d", &maxstep);
			break;
		case 'g':
			sscanf(optarg, "%lf", &tmp);
			gravity = (mtxfp)tmp;
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
