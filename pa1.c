
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <getopt.h>

#define MTXTYPEFLOAT
#include "matrix.h"

const mtxfp G = 0.0001;

struct particle {
	struct matrix *pos;
	mtxfp mass;
};

struct particle *initParticles(int pnum);

int runSim(int pnum, int maxstep, mtxfp timestep, mtxfp gravity)
{
	printf("Running with pnum: %d, maxstep: %d, timestep: %f, gravity: %f\n",
				 pnum, maxstep, timestep, gravity);
	struct particle *particles = initParticles(pnum);
	
	for(int i = 0; i < pnum; i++)
		mtxFree(particles[i].pos);
	free(particles);
	return 0;
}

struct particle *initParticles(int pnum)
{
	struct particle *particles = malloc(sizeof(struct particle[pnum]));
	assert(particles);
	for(int i = 0; i < pnum; i++) {
		particles[i].pos = mtxCreate(1, 3);
		mtxSet(particles[i].pos, 0, 0, ((mtxfp)rand()) / RAND_MAX);
		mtxSet(particles[i].pos, 0, 1, ((mtxfp)rand()) / RAND_MAX);
		mtxSet(particles[i].pos, 0, 2, ((mtxfp)rand()) / RAND_MAX);
    particles[i].mass = ((mtxfp)rand()) / RAND_MAX;
	}
	return particles;
}

int main(int argc, char **argv)
{
	struct matrix *mtx = mtxCreate(2, 3);
	mtxSet(mtx, 0, 0, (double)1.0);
	mtxSet(mtx, 0, 1, (double)10.0);
	mtxSet(mtx, 0, 2, (double)100.0);
	mtxSet(mtx, 1, 0, (double)5.0);
	mtxSet(mtx, 1, 1, (double)25.0);
	mtxSet(mtx, 1, 2, (double)125.0);
	mtxTranspose(mtx);
	mtxTranspose(mtx);
	mtxFree(mtx);
	float a, *b;
	b = &a;
	*b = 65.0;
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
