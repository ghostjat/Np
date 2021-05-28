#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <json-c/json.h>

#define _ME  2.7182818284590452354
#define _EPSILON  0.000000001
#define _PI  3.14159265358979323846
#define _TWOPI  6.28318530718
#define _SPACING 10
#define _FORMAT "%-10.3g"

#define max(a,b) \
({__typeof__ (a) _a = (a);\
__typeof__ (b) _b = (b); \
_a > _b ? _a : _b;})
#define min(a,b) \
({__typeof__ (a) _a = (a);\
__typeof__ (b) _b = (b); \
_a < _b ? _a : _b;})

typedef struct _matrix {
    int rows;
    int cols;
    double *data;
} matrix;

typedef struct _vector {
    int cols;
    double *data;
} vector;

typedef struct _lup {
    matrix *l;
    matrix *u;
    int *p;
    int *n;
} lup;

int isSquare(matrix *m);
int isSymmetric(matrix *m);
int checkShape(matrix *m1, matrix *m2);
int checkDimensions(matrix *m1, matrix *m2);