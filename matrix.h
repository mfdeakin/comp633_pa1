
/*
 * matrix.h
 * An implementation of basic matrix operations
 * By Michael Deakin
 * Still a work in progress
 */

#ifndef __MATRIX_H
#define __MATRIX_H

struct matrix;

enum matrixorder {
    MTXROWMAJORORDER,
    MTXCOLMAJORORDER
};

#ifdef MTXTYPEFLOAT
typedef float mtxfp;
#else
typedef double mtxfp;
#endif

struct matrix *mtxCreate(unsigned width, unsigned height);
struct matrix *mtxCreateI(unsigned size);
struct matrix *mtxCopy(struct matrix *mtx);
struct matrix *mtxFromArray(mtxfp array[], unsigned w, unsigned h);
void mtxFree(struct matrix *mtx);
void mtxSetOrder(struct matrix *mtx, enum matrixorder order);
struct matrix *mtxAdd(struct matrix *lhs, struct matrix *rhs);
struct matrix *mtxSub(struct matrix *lhs, struct matrix *rhs);
struct matrix *mtxNeg(struct matrix *mtx);
struct matrix *mtxMul(struct matrix *lhs, struct matrix *rhs);
struct matrix *mtxScalarMul(struct matrix *mtx, mtxfp scalar);
mtxfp mtxDeterminant(struct matrix *mtx);

mtxfp mtxGet(struct matrix *mtx, unsigned x, unsigned y);
int mtxSet(struct matrix *mtx, unsigned x, unsigned y, mtxfp val);
void mtxClear(struct matrix *mtx);

struct matrix *mtxTranslate(mtxfp x, mtxfp y, mtxfp z);
struct matrix *mtxScale(mtxfp x, mtxfp y, mtxfp z);
struct matrix *mtxRotation(mtxfp theta, mtxfp x, mtxfp y, mtxfp z);

#endif
