#include "numphp.h"

matrix *mat(int rows, int cols) {
    assert(rows > 0 && cols > 0);
    matrix *ar;
    ar = malloc(sizeof (matrix));
    ar->rows = rows;
    ar->cols = cols;
    ar->data = (double*) malloc((size_t) rows * cols * sizeof (double));
    return ar;
}

matrix *matNull(int rows, int cols) {
    matrix *ar = mat(rows, cols);
    return ar;
}

matrix *matZeros(int rows, int cols) {
    matrix *ar = mat(rows, cols);
    register int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = 0.0;
        }
    }
    return ar;
}

matrix *matOnes(int rows, int cols) {
    matrix *ar = mat(rows, cols);
    register int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = 1.0;
        }
    }
    return ar;
}

matrix *matFull(int rows, int cols, double val) {
    matrix *ar = mat(rows, cols);
    register int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = val;
        }
    }
    return ar;
}

matrix *matDiagonal(vector *v) {
    matrix *ar = mat(v->cols, v->cols);
    register int i;
    for (i = 0; i < v->cols; ++i) {
        ar->data[i * v->cols + i] = v->data[i];

    }
    return ar;
}

matrix *matPoisson(int rows, int cols, double lambda) {
    matrix *ar = mat(rows, cols);
    int max = RAND_MAX;
    double l = exp(-lambda);
    register int r, c;
    unsigned int k;
    double p;
    for (r = 0; r < rows; ++r) {
        for (c = 0; c < cols; ++c) {
            k = 0;
            p = 1.0;
            while (p > l) {
                k++;
                p = p * rand() / max;
            }
            ar->data[r * cols + c] = k - 1;
        }
    }
    return ar;
}

matrix *matMinimum(matrix *m1, matrix *m2) {
    if (checkShape(m1, m2)) {
        matrix *ar = mat(m1->rows, m2->cols);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * m2->cols + j] = min(m1->data[i * m2->cols + j], m2->data[i * m2->cols + j]);
            }
        }
        return ar;
    }
    exit(1);
}

matrix *matMaximum(matrix *m1, matrix *m2) {
    if (checkShape(m1, m2)) {
        matrix *ar = mat(m1->rows, m2->cols);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * m2->cols + j] = max(m1->data[i * m2->cols + j], m2->data[i * m2->cols + j]);
            }
        }
        return ar;
    }
    exit(1);
}

matrix *matRandom(int rows, int cols) {
    matrix *ar = mat(rows, cols);
    register int i, j;
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = 2 * (rand() / (double) RAND_MAX) - 1;
        }
    }
    return ar;
}

matrix *matUniform(int rows, int cols) {
    matrix *ar = mat(rows, cols);
    register int i, j;
    srand(time(NULL));
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = 2 * (rand() / (double) RAND_MAX) - 1;
        }
    }
    return ar;
}

matrix *matCopy(matrix *ar) {
    register int i, j;
    matrix *arCopy = mat(ar->rows, ar->cols);
    for (i = 0; i < ar->rows; ++i) {
        for (j = 0; j < ar->cols; ++j) {
            arCopy->data[i * ar->cols + j] = ar->data[i * ar->cols + j];
        }
    }
    return arCopy;
}

matrix *matIdentity(int rows, int cols) {
    register int i, j;
    matrix *ar = mat(rows, cols);
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j) {
            ar->data[i * cols + j] = i == j ? 1 : 0;
        }
    }
    return ar;
}

matrix *matMultiply(matrix* m1, matrix* m2) {
    assert(m1->rows == m2->rows && m1->cols == m2->cols);
    matrix *ar = mat(m1->rows, m2->cols);
    register int i, j;
    for (i = 0; i < m1->rows; ++i) {
        for (j = 0; j < m2->cols; ++j) {
            ar->data[i * m2->cols + j] = m1->data[i * m2->cols + j] * m2->data[i * m2->cols + j];
        }
    }
    return ar;
}

void dignoalInterChange(matrix *m) {
    register int i, j;
    double tmp;
    for (i = 0; i < m->rows; ++i) {
        for (j = 0; j < m->cols; ++j) {
            tmp = m->data[i * m->cols - j];
            m->data[i * m->cols - j] = m->data[i * m->cols + j];
            m->data[i * m->cols + j] = tmp;
        }
    }
}

double matTrace(matrix *m) {
    if (isSquare(m)) {
        register int i, j;
        double trace = 0.0;
        for (i = 0; i < m->rows; ++i) {
            for (j = 0; j < m->cols; ++j) {
                if (i == j) {
                    trace = trace + m->data[i * m->cols + i];
                }
            }
        }
        return trace;
    }
    exit(1);
}

double matCell(matrix *m, int row, int col) {
    assert(col < m->cols && row < m->rows);
    return m->data[row * m->cols + col];
}

void matPrint(matrix* mat) {
    register int i, j;
    for (i = 0; i < mat->rows; ++i) {
        printf("%s", "[");
        for (j = 0; j < mat->cols; ++j) {
            printf("%lf \t", mat->data[i * mat->cols + j]);
        }
        printf("%s\n", "]");
    }
}


//--------------Vector functions-------------------

vector *vec(int cols) {
    vector *v;
    v = malloc(sizeof (vector));
    v->cols = cols;
    v->data = (double*) malloc((size_t) (cols * sizeof (double)));
    if (!v->data) {
        printf("vector: allocation error!\n");
        exit(1);
    } else {
        return v;
    }
}

vector *flatten(matrix *m) {
    vector *v = vec(m->rows * m->cols);
    register int i, j;
    for (i = 0; i < m->rows; ++i) {
        for (j = 0; j < m->cols; ++j) {
            v->data[i * m->cols + j] = m->data[i * m->cols + j];
        }
    }
    return v;
}

vector *matRowAsVector(matrix *m, int i) {
    vector *v = vec(m->cols);
    register int j;
    for (j = 0; j < m->cols; ++j) {
        v->data[j] = m->data[i * m->cols + j];
    }
    return v;
}

vector *matColAsVector(matrix *m, int j) {
    vector *v = vec(m->rows);
    register int i;
    for (i = 0; i < m->rows; ++i) {
        v->data[i] = m->data[i * m->rows + j];
    }
    return v;
}

vector *matDiagonalAsVector(matrix *m) {
    if (isSquare(m)) {
        vector *v = vec(m->rows);
        register int i;
        for (i = 0; i < m->rows; ++i) {
            v->data[i] = m->data[i * m->rows + i];
        }
        return v;
    } else {
        printf("%s\n", "Error::Can not trace of a none square matrix!");
        exit(1);
    }
}

//-----------------------End Vector----------------


//---------Matrix Opreations functions-------------------

matrix *matDgemm(matrix* m1, matrix* m2) {
    matrix *m3 = mat(m1->rows, m2->cols);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m1->rows, m2->rows, m1->cols, 1.0, m1->data, m2->cols, m2->data, m2->rows, 0.0, m3->data, m3->rows);
    return m3;
}

matrix *matInverse(matrix *m) {
    if (isSquare(m)) {
        matrix * ar = matCopy(m);
        int *ipiv = (int*) malloc((size_t) (m->cols * sizeof (int)));
        int lp = LAPACKE_dgetrf(101, ar->rows, ar->cols, ar->data, ar->rows, ipiv);
        if (lp != 0) {
            printf("%s\n", "Error::LAPACK numphp.c");
            exit(1);
        }

        lp = LAPACKE_dgetri(101, ar->rows, ar->data, ar->rows, ipiv);
        if (lp != 0) {
            printf("%s\n", "Error::LAPACK numphp.c");
            exit(1);
        }
        free(ipiv);
        return ar;
    } else {
        printf("%s\n", "Error::given matris is not a Square matrix!");
        exit(1);
    }
}

matrix *matRef(matrix *m) {
    int *ipiv = (int*) malloc((size_t) (m->cols * sizeof (int)));
    matrix *ar = matCopy(m);
    int lp = LAPACKE_dgetrf(101, ar->rows, ar->cols, ar->data, ar->rows, ipiv);
    if (lp != 0) {
        printf("%s\n", "Error::LAPACK  numphp.c");
        exit(1);
    }
    free(ipiv);
    return ar;

}

matrix *matCholesky(matrix *m) {
    if (isSymmetric(m)) {
        matrix *ar = matCopy(m);
        int lp = LAPACKE_dpotrf(101, 'L', ar->cols, ar->data, ar->cols);
        if (lp != 0) {
            printf("%s\n", "Error::LAPACK numphp.c");
            exit(1);
        }
        register int i, j;
        for (i = 0; i < ar->cols; ++i) {
            for (j = i + 0; j < ar->cols; ++j) {
                ar->data[i * ar->cols + j] = 0.0;
            }
        }
        return ar;
    }
    printf("%s\n", "Error::given matrix is not a Symmetric!");
    exit(1);
}

vector *matArgMax(matrix *m) {
    register int i;
    vector *v = vec(m->rows);
    for (i = 0; i < m->rows; ++i) {
        vector *vr = matRowAsVector(m, i);
        v->data[i] = cblas_idamax(vr->cols, vr->data, 1);
    }
    return v;
}

vector *matArgMin(matrix *m) {
    register int i;
    vector *v = vec(m->rows);
    for (i = 0; i < m->rows; ++i) {
        vector *vr = matRowAsVector(m, i);
        v->data[i] = cblas_idamin(vr->cols, vr->data, 1);
    }
    return v;
}

matrix *matAdd(matrix *m1, matrix *m2) {
    if (checkShape(m1, m2)) {
        matrix *ar = mat(m1->rows, m2->cols);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * m2->cols + j] = m1->data[i * m2->cols + j] + m2->data[i * m2->cols + j];
            }
        }
        return ar;
    } else {
        exit(1);
    }
}

matrix *matTranspose(matrix *m) {
    matrix *at = mat(m->cols, m->rows);
    register int i, j;
    for (i = 0; i < m->rows; ++i) {
        for (j = 0; j < m->cols; ++j) {
            at->data[i * m->cols + j] = m->data[j * m->cols + i];
        }
    }
    return at;
}

matrix *matJoinLeft(matrix *m1, matrix *m2) {
    if (m1->rows == m2->rows) {
        int cols = m1->cols + m2->cols;
        matrix *ar = mat(m1->rows, cols);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m1->cols; ++j) {
                ar->data[i * cols + j] = m1->data[i * m1->cols + j];
            }
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * cols + (m1->cols + j)] = m2->data[i * m2->cols + j];
            }
        }
        return ar;
    }
    else {
        printf("%s\n", "Error::Invalid size!");
        exit(1);
    }
}

matrix *matJoinRight(matrix *m1, matrix *m2) {
    if (m1->rows == m2->rows) {
        int col = m1->cols + m2->cols;
        matrix *ar = mat(m1->rows, col);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m1->cols; ++j) {
                ar->data[i * m1->cols + j] = m1->data[i * m1->cols + j];
            }
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * col + (m2->cols + j)] = m2->data[i * m2->cols + j];
            }
        }
        return ar;
    } else {
        printf("%s\n", "Error::Invalid size!");
        exit(1);
    }

}

matrix *matJoinAbove(matrix *m1, matrix *m2) {
    if (m1->cols == m2->cols) {
        int row = m1->rows + m2->rows;
        matrix *ar = mat(row, m1->cols);
        register int i, j;
        for (i = 0; i < m2->rows; ++i) {
            for (j = 0; j < m2->cols; ++j) {
                ar->data[i * m2->cols + j] = m2->data[i * m2->cols + j];
            }
            for (j = 0; j < m1->cols; ++j) {
                ar->data[(i + m1->rows) * m1->cols + j] = m1->data[i * m1->cols + j];
            }
        }
        return ar;
    } else {
        printf("%s\n", "Error::Invalid size!");
        exit(1);
    }

}

matrix *matJoinBelow(matrix *m1, matrix *m2) {
    if (m1->cols == m2->cols) {
        int row = m1->rows + m2->rows;
        matrix *ar = mat(row, m1->cols);
        register int i, j;
        for (i = 0; i < m1->rows; ++i) {
            for (j = 0; j < m1->cols; ++j) {
                ar->data[i * m1->cols + j] = m1->data[i * m1->cols + j];
            }
            for (j = 0; j < m2->cols; ++j) {
                ar->data[(i + m2->rows) * m2->cols + j] = m2->data[i * m2->cols + j];
            }
        }
        return ar;
    } else {
        printf("%s\n", "Error::Invalid size!");
        exit(1);
    }

}

//---------------check type functions------------

int checkShape(matrix *m1, matrix *m2) {
    if (m1->rows == m2->rows && m1->cols == m2->cols) {
        return 1;
    }

    printf("%s\n", "Error::mismatch shape of given matrix!");
    return 0;
}

int checkDimensions(matrix *m1, matrix *m2) {
    if (m1->rows == m2->cols) {
        return 1;
    }

    printf("%s\n", "Error::mismatch dimensions of given matrix!");
    return 0;
}

int isSquare(matrix *m1) {
    if (m1->rows == m1->cols) {
        return 1;
    }
    printf("%s\n", "Error::given matrix is not squared!");
    return 0;
}

int isSymmetric(matrix *m) {
    if (!isSquare(m)) {
        exit(1);
    }
    matrix *ar = matTranspose(m);
    register int i, j;
    for (i = 0; i < m->rows; ++i) {
        for (j = 0; j < m->cols; ++j) {
            if (ar->data[i * m->cols + j] != m->data[i * m->cols + j]) {
                free(ar);
                return 0;
            }
        }
    }
    free(ar);
    return 1;
}

int isZero(double d) {
    return fabs(d) < _EPSILON;
}

//---------------------End Check Type--------------------------

//-----------------------general functions----------------------

void freeMatrix(matrix* m) {
    if (m != NULL) {
        if (m->data != NULL) {
            free(m->data);
            m->data = NULL;
        }
        free(m);
        m = NULL;
    }
    return;
}