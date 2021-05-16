#define FFI_SCOPE "lapack"
#define FFI_LIB "liblapacke.so"

int LAPACKE_sgetrf( int matrix_layout, int m, int n,
                           float* a, int lda, int* ipiv );
int LAPACKE_sgetri( int matrix_layout, int n, float* a,
                           int lda, const int* ipiv );
int LAPACKE_sgesdd( int matrix_layout, char jobz, int m,
                           int n, float* a, int lda, float* s,
                           float* u, int ldu, float* vt,
                           int ldvt );
int LAPACKE_sgeev( int matrix_layout, char jobvl, char jobvr,
                          int n, float* a, int lda, float* wr,
                          float* wi, float* vl, int ldvl, float* vr,
                          int ldvr );
int LAPACKE_ssyev( int matrix_layout, char jobz, char uplo, int n,
                          float* a, int lda, float* w );
int LAPACKE_sgels( int matrix_order, char trans, int m,
                          int n, int nrhs, float* a,
                          int lda, float* b, int ldb );
int LAPACKE_sbdsdc( int matrix_layout, char uplo, char compq,
                           int n, float* d, float* e, float* u,
                           int ldu, float* vt, int ldvt, float* q,
                           int* iq );
int LAPACKE_sbdsvdx( int matrix_layout, char uplo, char jobz, char range,
                           int n, float* d, float* e,
                           float vl, float vu,
                           int il, int iu, int* ns,
                           float* s, float* z, int ldz,
                           int* superb );
int LAPACKE_sdisna( char job, int m, int n, const float* d,
                           float* sep );
int LAPACKE_sgbbrd( int matrix_layout, char vect, int m,
                           int n, int ncc, int kl,
                           int ku, float* ab, int ldab, float* d,
                           float* e, float* q, int ldq, float* pt,
                           int ldpt, float* c, int ldc );
int LAPACKE_sgbcon( int matrix_layout, char norm, int n,
                           int kl, int ku, const float* ab,
                           int ldab, const int* ipiv, float anorm,
                           float* rcond );
int LAPACKE_sgbequ( int matrix_layout, int m, int n,
                           int kl, int ku, const float* ab,
                           int ldab, float* r, float* c, float* rowcnd,
                           float* colcnd, float* amax );
int LAPACKE_sgbequb( int matrix_layout, int m, int n,
                            int kl, int ku, const float* ab,
                            int ldab, float* r, float* c, float* rowcnd,
                            float* colcnd, float* amax );
int LAPACKE_sgbrfs( int matrix_layout, char trans, int n,
                           int kl, int ku, int nrhs,
                           const float* ab, int ldab, const float* afb,
                           int ldafb, const int* ipiv,
                           const float* b, int ldb, float* x,
                           int ldx, float* ferr, float* berr );
int LAPACKE_sgbsv( int matrix_layout, int n, int kl,
                          int ku, int nrhs, float* ab,
                          int ldab, int* ipiv, float* b,
                          int ldb );
int LAPACKE_sgbsvx( int matrix_layout, char fact, char trans,
                           int n, int kl, int ku,
                           int nrhs, float* ab, int ldab,
                           float* afb, int ldafb, int* ipiv,
                           char* equed, float* r, float* c, float* b,
                           int ldb, float* x, int ldx,
                           float* rcond, float* ferr, float* berr,
                           float* rpivot );
int LAPACKE_sgbtrf( int matrix_layout, int m, int n,
                           int kl, int ku, float* ab,
                           int ldab, int* ipiv );
int LAPACKE_sgbtrs( int matrix_layout, char trans, int n,
                           int kl, int ku, int nrhs,
                           const float* ab, int ldab,
                           const int* ipiv, float* b, int ldb );
int LAPACKE_sgebak( int matrix_layout, char job, char side, int n,
                           int ilo, int ihi, const float* scale,
                           int m, float* v, int ldv );
int LAPACKE_sgebal( int matrix_layout, char job, int n, float* a,
                           int lda, int* ilo, int* ihi,
                           float* scale );
int LAPACKE_sgebrd( int matrix_layout, int m, int n,
                           float* a, int lda, float* d, float* e,
                           float* tauq, float* taup );
int LAPACKE_sgecon( int matrix_layout, char norm, int n,
                           const float* a, int lda, float anorm,
                           float* rcond );
int LAPACKE_sgeequ( int matrix_layout, int m, int n,
                           const float* a, int lda, float* r, float* c,
                           float* rowcnd, float* colcnd, float* amax );
int LAPACKE_sgeequb( int matrix_layout, int m, int n,
                            const float* a, int lda, float* r, float* c,
                            float* rowcnd, float* colcnd, float* amax );
int LAPACKE_sgelq2( int matrix_layout, int m, int n,
                           float* a, int lda, float* tau );
int LAPACKE_spotrf( int matrix_order, char uplo, int n, float* a,
                           int lda );
int LAPACKE_sgelsy( int matrix_order, int m, int n, int nrhs, float* a, int lda, float* b,
                           int ldb, int* jpvt, float rcond, int* rank );
void LAPACKE_sgelsd( int m, int n, int nrhs, float* a,
                    int lda, float* b, int ldb, float* s,
                    float* rcond, int rank, float* work,
                    int lwork, int iwork, int info );
float LAPACKE_slange(int matrix_layout, char norm, int m,
                        int n, const float* a, int lda);
