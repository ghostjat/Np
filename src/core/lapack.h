#define FFI_SCOPE "lapack"
#define FFI_LIB "liblapacke.so"


int LAPACKE_dgetrf( int matrix_layout, int m, int n,
                           double* a, int lda, int* ipiv );

int LAPACKE_dgetri( int matrix_layout, int n, double* a,
                           int lda, const int* ipiv );

int LAPACKE_dgesdd( int matrix_layout, char jobz, int m,
                           int n, double* a, int lda, double* s,
                           double* u, int ldu, double* vt,
                           int ldvt );

int LAPACKE_dgeev( int matrix_layout, char jobvl, char jobvr,
                          int n, double* a, int lda, double* wr,
                          double* wi, double* vl, int ldvl, double* vr,
                          int ldvr );

int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, int n,
                          double* a, int lda, double* w );

int LAPACKE_dgels( int matrix_order, char trans, int m,
                          int n, int nrhs, double* a,
                          int lda, double* b, int ldb );
int LAPACKE_dgebal( int matrix_layout, char job, int n, double* a,
                           int lda, int* ilo, int* ihi,double* scale );

int LAPACKE_dpotrf( int matrix_layout, char uplo, int n, double* a,
                           int lda );

float LAPACKE_dlange(int matrix_layout, char norm, int m,
                        int n, const double* a, int lda);


int LAPACKE_dlasrt( char id, int n, double* d );

