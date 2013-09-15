
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

struct matrix {
	/* For easy switching between matrix implementations */
	mtxfp * (*matrixPos)(struct matrix *mtx, unsigned x, unsigned y);
	unsigned width, height;
	mtxfp mtx[];
};

#ifdef MTXDEBUG
#define debug(str) fprintf(stderr, str)

void mtxPrint(struct matrix *m)
{
	m->matrixPos(m, 0, 0);
	for(unsigned i = 0, j = 0; i < m->height;
			j = (j + 1) % m->width, i += j == 0) {
		if(j == 0)
			printf("\n");
		printf("%4f ", *(m->matrixPos(m, j, i)));
	}
	printf("\n");
}
#else
int debug(char *str) {
	return 0;
}

#define mtxPrint(m)
#endif

#if MTXMAJORORDER == MTXMAJORROW
#define MTXMAJORFUNC matrixPosRow
#define MTXSTRIDE width
#define MTXRUN height
#else
#define MTXMAJORFUNC matrixPosCol
#define MTXSTRIDE height
#define MTXRUN width
#endif

/* Calculates the position in the matrix array of the value */
mtxfp *matrixPosCol(struct matrix *mtx, unsigned x, unsigned y)
{
	assert(mtx && x < mtx->width && y < mtx->height);
	/* For column major order: */
	return mtx->mtx + (mtx->height * x + y);
}

mtxfp *matrixPosRow(struct matrix *mtx, unsigned x, unsigned y)
{
	assert(mtx && x < mtx->width && y < mtx->height);
	/* For row major order: */
	return mtx->mtx + (mtx->width * y + x);
}

struct matrix *mtxCreate(unsigned width, unsigned height)
{
	if(!width || !height) {
		debug("mtxCreate: Width and height are zero\n");
		return NULL;
	}
	struct matrix *mtx = malloc(sizeof(struct matrix) + sizeof(mtxfp[width][height]));
	assert(mtx);
	memset(mtx->mtx, 0.0, sizeof(mtxfp[width][height]));
	mtx->width = width;
	mtx->height = height;
	mtx->matrixPos = MTXMAJORFUNC;
	return mtx;
}

struct matrix *mtxCreateI(unsigned size)
{
	if(!size) {
		debug("mtxCreateI: size is zero\n");
		return NULL;
	}
	struct matrix *mtx = malloc(sizeof(struct matrix) + sizeof(mtxfp[size][size]));
	assert(mtx);
	memset(mtx->mtx, 0.0f, sizeof(mtxfp[size][size]));
	mtx->width = size;
	mtx->height = size;
	for(int i = 0; i < mtx->MTXRUN * mtx->MTXSTRIDE;
			i += mtx->MTXSTRIDE + 1)
		mtx->mtx[i] = 1.0f;
	mtx->matrixPos = MTXMAJORFUNC;
	return mtx;
}

struct matrix *mtxCopy(struct matrix *mtx)
{
	assert(mtx);
	struct matrix *copy = mtxCreate(mtx->width, mtx->height);
	int i;
	for(i = 0; i < (mtx->width * mtx->height); i++)
		copy->mtx[i] = mtx->mtx[i];
	return copy;
}

struct matrix *mtxFromArray(mtxfp array[], unsigned w, unsigned h)
{
	assert(array);
	struct matrix *mtx = mtxCreate(w, h);
	unsigned i, x, y;
	for(i = 0, x = 0, y = 0;
			i < w * h;
			i++, x++) {
		if(x == mtx->width) {
			x = 0;
			y++;
		}
		*mtx->matrixPos(mtx, x, y) = array[i];
	}
	return mtx;
}

void mtxFree(struct matrix *mtx)
{
	assert(mtx);
	free(mtx);
}

void mtxSetOrder(struct matrix *mtx, enum matrixorder order)
{
	assert(mtx &&
				 (order == MTXROWMAJORORDER ||
					order == MTXCOLMAJORORDER));
	if(order == MTXROWMAJORORDER)
		mtx->matrixPos = matrixPosRow;
	else
		mtx->matrixPos = matrixPosCol;
}

void mtxTranspose(struct matrix *mtx)
{
	assert(mtx);
	if(mtx->matrixPos == matrixPosRow) {
		mtx->matrixPos = matrixPosCol;
		unsigned tmp = mtx->width;
		mtx->width = mtx->height;
		mtx->height = tmp;
	}
	else {
		mtx->matrixPos = matrixPosRow;
		unsigned tmp = mtx->width;
		mtx->width = mtx->height;
		mtx->height = tmp;
	}
}

struct matrix *mtxAdd(struct matrix *lhs, struct matrix *rhs)
{
	assert(lhs && rhs);
	if(lhs->width != rhs->width || lhs->height != rhs->height) {
		debug("mtxAdd: lhs and rhs are not valid for adding\n");
		return NULL;
	}
	struct matrix *mtx = mtxCreate(lhs->width, lhs->height);
	unsigned x = 0, y = 0;
	for(unsigned i = 0; i < (lhs->width * lhs->height); i++, x++) {
		if(x == lhs->width) {
			x = 0;
			y++;
		}
		*mtx->matrixPos(mtx, x, y) = *lhs->matrixPos(lhs, x, y) + *rhs->matrixPos(rhs, x, y);
	}
	return mtx;
}

struct matrix *mtxSub(struct matrix *lhs, struct matrix *rhs)
{
	assert(lhs && rhs);
	if(lhs->width != rhs->width || lhs->height != rhs->height) {
		debug("mtxSub: lhs and rhs are not valid for subtracting\n");
		return NULL;
	}
	struct matrix *mtx = mtxCreate(lhs->width, lhs->height);
	unsigned x = 0, y = 0;
	for(unsigned i = 0; i < (lhs->width * lhs->height); i++, x++) {
		if(x == lhs->width) {
			x = 0;
			y++;
		}
		*mtx->matrixPos(mtx, x, y) = *lhs->matrixPos(lhs, x, y) -
			*rhs->matrixPos(rhs, x, y);
	}
	return mtx;
}

struct matrix *mtxNeg(struct matrix *mtx)
{
	assert(mtx);
	struct matrix *ret = mtxCreate(mtx->width, mtx->height);
	for(unsigned i = 0; i < (mtx->width * mtx->height); i++)
		ret->mtx[i] = -mtx->mtx[i];
	return ret;
}

struct matrix *mtxMul(struct matrix *lhs, struct matrix *rhs)
{
	assert(lhs && rhs);
	if(lhs->width != rhs->height) {
		debug("mtxMul: lhs and rhs cannot be multiplied\n");
		return NULL;
	}
	struct matrix *mtx = mtxCreate(rhs->width, lhs->height);
	unsigned x = 0, y = 0;
	for(unsigned i = 0; i < (mtx->width * mtx->height); i++, x++) {
		if(x == mtx->width) {
			x = 0;
			y++;
		}
		for(unsigned xypos = 0; xypos < lhs->width; xypos++) {
			*mtx->matrixPos(mtx, x, y) += *lhs->matrixPos(lhs, xypos, y) *
				*rhs->matrixPos(rhs, x, xypos);
		}
	}
	return mtx;
}

struct matrix *mtxScalarMul(struct matrix *mtx, mtxfp scalar)
{
	assert(mtx);
	struct matrix *ret = mtxCreate(mtx->width, mtx->height);
	for(unsigned i = 0; i < (mtx->width * mtx->height); i++)
		ret->mtx[i] = mtx->mtx[i] * scalar;
	return ret;
}

mtxfp mtxDeterminant(struct matrix *mtx)
{
	assert(mtx);
	if(mtx->width != mtx->height) {
		debug("mtxDeterminate: Not a square matrix\n");
		/* Return NaN */
		return 0.0 / 0.0;
	}
	if(mtx->width == 1) {
		return mtx->mtx[0];
	}
	if(mtx->width == 2) {
		return *mtx->matrixPos(mtx, 0, 0) * *mtx->matrixPos(mtx, 1, 1) -
			*mtx->matrixPos(mtx, 1, 0) * *mtx->matrixPos(mtx, 0, 1);
	}
	struct matrix *buffer = mtxCreate(mtx->width - 1, mtx->height - 1);
	mtxfp determinant = 0;
	signed sign = 1;
	for(unsigned i = 0; i < mtx->width; i++) {
		unsigned outOffset = 0;
		for(unsigned pos = 0; pos < mtx->width; pos++) {
			if(pos == i)
				continue;
			mtxfp *inOffset = mtx->matrixPos(mtx, pos, 1);
			for(unsigned j = 1;
					j < mtx->width;
					j++, inOffset++, outOffset++) {
				buffer->mtx[outOffset] = *inOffset;
			}
		}
		determinant += sign * mtx->mtx[i] * mtxDeterminant(buffer);
		sign = -sign;
	}
	return determinant;
}

mtxfp mtxGet(struct matrix *mtx, unsigned x, unsigned y)
{
	if(x >= mtx->width || y >= mtx->height) {
		debug("mtxGet: Invalid index\n");
		return 0.0;
	}

	return *mtx->matrixPos(mtx, x, y);
}

int mtxSet(struct matrix *mtx, unsigned x, unsigned y, mtxfp val)
{
	if(x >= mtx->width || y >= mtx->height) {
		debug("mtxSet: Invalid index\n");
		return -1;
	}
	*mtx->matrixPos(mtx, x, y) = val;
	return 0;
}

void mtxClear(struct matrix *mtx)
{
	for(int i = 0; i < mtx->width * mtx->height; i++)
		mtx->mtx[i] = 0;
}

struct matrix *mtxRotation(mtxfp theta, mtxfp x, mtxfp y, mtxfp z)
{
	mtxfp crot = cos(theta), srot = sin(theta);
	mtxfp values[16] = {
		crot + x * x * (1 - crot), x * y * (1 - crot) - z * srot, x * z * (1 - crot) + y * srot, 0.0f,
		y * x * (1 - crot) + z * srot, crot + y * y * (1 - crot), y * z * (1 - crot) - x * srot, 0.0f,
		z * x * (1 - crot) - y * srot, z * y * (1 - crot) + x * srot, crot + z * z * (1 - crot), 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	return mtxFromArray(values, 4, 4);
}

struct matrix *mtxTranslate(mtxfp x, mtxfp y, mtxfp z)
{
	mtxfp values[16] = {
		1.0f, 0.0f, 0.0f, x,
		0.0f, 1.0f, 0.0f, y,
		0.0f, 0.0f, 1.0f, z,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	return mtxFromArray(values, 4, 4);
}

struct matrix *mtxScale(mtxfp x, mtxfp y, mtxfp z)
{
	mtxfp values[16] = {
		x, 0.0f, 0.0f, 0.0f,
		0.0f, y, 0.0f, 0.0f,
		0.0f, 0.0f, z, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	return mtxFromArray(values, 4, 4);
}
