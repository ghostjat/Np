<?php

namespace numphp\core;

/**
 * php interface for OpenBLAS
 * 
 * @package NumPhp\Blas
 * @category Scientific Computing
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class blas {

    const CblasLeft = 141, CblasRight = 142;
    const CblasUpper = 121, CblasLower = 122;
    const CblasNonUnit = 131, CblasUnit = 132;
    const CblasRowMajor = 101, CblasColMajor = 102;
    const CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113;

    public static $ffi_blas = null;

    
    public static function init() {
        if (is_null(self::$ffi_blas)) {
            self::$ffi_blas = \FFI::load(__DIR__ . '/blas.h');
        }
        return self::$ffi_blas;
    }
    
    public static function setNumThreads(int $num_threads) {
        self::init();
        self::$ffi_blas->openblas_set_num_threads($num_threads);
    }
    
    public static function getNumThreads() :int {
        self::init();
        return self::$ffi_blas->openblas_get_num_threads();
    }
    
    public static function getNumPorcs() :int {
        self::init();
        return self::$ffi_blas->openblas_get_num_procs();
    }
    
    public static function getConfig() {
        self::init();
        return self::$ffi_blas->openblas_get_config();
    }
    
    public static function getCoreName() {
        self::init();
        return self::$ffi_blas->openblas_get_corename();
    }
    
    public static function getNumParallel() :int {
        self::init();
        return self::$ffi_blas->openblas_get_parallel();
    }

    /**
     * Product of general matrix and general matrix
     *    C := alpha * AB + beta * C
     * @dtype Float
     * @param \numphp\matrix $m1
     * @param \numphp\matrix $m2
     * @param \numphp\matrix $mr
     * @return \FFI\CData
     */
    public static function sgemm(\numphp\matrix $m1, \numphp\matrix $m2, \numphp\matrix $mr,int $trans1= self::CblasNoTrans,int $trans2= self::CblasNoTrans) {
        self::init();
        self::$ffi_blas->cblas_sgemm(self::CblasRowMajor, $trans1, $trans2, $m1->row, $m2->col, $m1->col, 1.0, $m1->data, $m1->col, $m2->data, $m2->col, 0.0, $mr->data, $mr->col);
    }

    /**
     *  Product of symmetric matrix and general matrix
     *    C := alpha * AB + beta * C
     *            or
     *    C := alpha * BA + beta * C
     * 
     * @param \numphp\matrix $m1
     * @param \numphp\matrix $m2
     * @param \numphp\matrix $mr
     */
     
    public static function ssymm(\numphp\matrix $m1, \numphp\matrix $m2, \numphp\matrix $mr) {
        self::init();
        self::$ffi_blas->cblas_ssymm(self::CblasRowMajor, self::CblasLeft, self::CblasUpper, $m1->row, 
                $m2->col,1.0, $m1->data, $m1->row, $m2->data, $m2->row, 0.0, $mr->data, $mr->row);
    }
    
    /**
     * Update rank n of symmetric matrix
     *    C := alpha * A A^T + beta * C
     *            or
     *    C := alpha * A^T A + beta * C
     * 
     * @param \numphp\matrix $m1
     * @param \numphp\matrix $m2
     */
    public static function ssyrk(\numphp\matrix $m1, \numphp\matrix $m2) {
        self::init();
        self::$ffi_blas->cblas_ssyrk(self::CblasRowMajor, self::CblasUpper,
        self::CblasNoTrans, $m1->row, $m2->col, 1.0, $m1->data, $m1->row, 0.0, $m2->data, $m2->row);
    }
    
    /**
     * Update rank 2k of symmetric matrix
     *    C := alpha * A B^T + alpha B A^T + beta * C
     *            or
     *    C := alpha * B^T A + alpha A^T B + beta * C
     * 
     * @param \numphp\matrix $m1
     * @param \numphp\matrix $m2
     * @param \numphp\matrix $mr
     */
    public static function ssyr2k(\numphp\matrix $m1, \numphp\matrix $m2, \numphp\matrix $mr) {
        self::init();
        self::$ffi_blas->cblas_ssyr2k(self::CblasRowMajor, self::CblasLower, self::CblasNoTrans, $m1->col, $m2->row, 1.0, $m1->data, $m1->row,
                  $m2->data, $m2->row, 0.0, $mr->data, $mr->row);
    }


    /**
     * @static
     * @dtype Double
     * @param \numphp\matrix $m1
     * @param \numphp\matrix $m2
     * @param \numphp\matrix $mr
     * @return \FFI\CData
     */
    public static function dgemm(\numphp\matrix $m1, \numphp\matrix $m2, \numphp\matrix $mr) {
        self::init();
        self::$ffi_blas->cblas_dgemm(self::CblasRowMajor, self::CblasNoTrans, self::CblasNoTrans, $m1->row,
                $m2->col, $m1->col, 1.0, $m1->data, $m1->col, $m2->data, $m2->row, 0.0, $mr->data, $mr->row);
    }

    /**
     * @static
     * @dtype Double
     * @param \numphp\matrix $m
     * @param \numphp\vector $v
     * @param \numphp\vector $mvr
     * @return \FFI\CData
     */
    public static function dgemv(\numphp\matrix $m, \numphp\vector $v, \numphp\vector $mvr) {
        self::init();
        self::$ffi_blas->cblas_dgemv(self::CblasRowMajor, self::CblasNoTrans, $m->row, $m->col,
                1.0, $m->data, $m->row, $v->data, 1, 1.0, $mvr->data, 1);
    }

    /**
     *  Compute the product of a general matrix and a vector.
     *    y := alpha * Ax + beta * y
     * @dtype Float
     * @param \numphp\matrix $m
     * @param \numphp\vector $v
     * @param \numphp\vector $mvr
     * @return \FFI\CData
     */
    public static function sgemv(\numphp\matrix $m, \numphp\vector $v, \numphp\vector $mvr) {
        self::init();
        return self::$ffi_blas->cblas_sgemv(self::CblasRowMajor, self::CblasNoTrans, $m->row, $m->col,
                1.0, $m->data, $m->row, $v->data, 1, 1.0, $mvr->data, 1);
    }
    
    /**
     * Compute the product of a general matrix and a vector stored in band format.
     *    y := alpha * Ax + beta * y
     * @param int $KL   Number of elements in the lower left part
     * @param int $KU   Number of elements in the upper right part
     * @param float $alpha  Coefficient of scalar multiple of vector
     * @param float $beta   Coefficient of scalar multiple of  mvr
     * @param \numphp\matrix $matrix
     * @param \numphp\vector $vector
     * @param \numphp\vector $mvr
     */
    public static function sgbmv(int $KL, int $KU, float $alpha, float $beta, \numphp\matrix $matrix, \numphp\vector $vector, \numphp\vector $mvr) {
        self::init();
        return self::$ffi_blas->cblas_sgbmv(self::CblasRowMajor,
                 self::CblasNoTrans, $matrix->row, $matrix->col,
                 $KL, $KU, $alpha,
                 $matrix->data,$matrix->row, $vector->data,
                 1, $beta, $mvr->data, 1);

    }
    
    /**
     * Compute the product of a column vector and a row vector. (Real number)
     *    A := alpha * x y^t + A
     * @param \numphp\vector $v1
     * @param \numphp\vector $v2
     * @param \numphp\matrix $m
     * @return void
     */
    public static function sger(\numphp\vector $v1, \numphp\vector $v2, \numphp\matrix $m) {
        self::init();
        return self::$ffi_blas->cblas_sger(self::CblasRowMajor, $m->row, $m->col,
                1.0, $v1->data, 1, $v2->data, 1, $m->data, $m->row);
    }
    

    /**
     * @static
     * @dtype Double
     * @param \numphp\vector $v1
     * @param \numphp\vector $v2
     */
    public static function ddot(\numphp\vector $v1, \numphp\vector $v2) {
        self::init();
        return self::$ffi_blas->cblas_ddot($v1->col, $v1->data, 1,$v2->data, 1);
    }
    
    /**
     * Calculate the inner product value between vectors.
     * @dtype Float
     * @param \numphp\vector $v1
     * @param \numphp\vector $v2
     * @return float
     */
    public static function sdot(\numphp\vector $v1, \numphp\vector $v2) {
        self::init();
        return self::$ffi_blas->cblas_sdot($v1->col, $v1->data, 1,$v2->data, 1);
    }
    
    /**
     * @satatic
     * @param \numphp\vector $v
     * 
     */
    public static function dmax(\numphp\vector $v) {
        self::init();
        return self::$ffi_blas->cblas_isamax($v->col, $v->data, 1);
    }
    
    /**
     * Calculates the index of the element with the largest absolute value in the vector.
     *  Note that this subscript starts from 1. If 0 is returned, n is invalid.
     *    ret := arg max |X(i)|
     * @param \numphp\vector $v
     * @return int index of the element(Note that start from 0 according to cblas_)
     */
    public static function smax(\numphp\vector $v) {
        self::init();
        return self::$ffi_blas->cblas_isamax($v->col, $v->data, 1);
    }

    /**
     * @static
     * @param \numphp\vector $v
     * 
     */
    public static function dmin(\numphp\vector $v) {
        self::init();
        return self::$ffi_blas->cblas_idamin($v->col, $v->data, 1);
    }
    
    /**
     *  Calculates the index of the element with the smallest absolute value in the vector.
     *  Note that this subscript starts from 1. If 0 is returned, n is invalid.
     *    ret := arg min |X(i)|
     * @param \numphp\vector $v
     * @return int
     */
    public static function smin(\numphp\vector $v) {
        self::init();
        return self::$ffi_blas->cblas_isamin($v->col, $v->data, 1);
    }
    
    public static function dswap() {
        self::init();
        
    }
    /**
     * Exchange the contents of the vector.
     *    X := Y
     *    Y := X
     * @param \numphp\vector $v1
     * @param \numphp\vector $v2
     */
    public static function sswap(\numphp\vector $v1, \numphp\vector $v2) {
        self::init();
        self::$ffi_blas->cblas_sswap($v1->col, $v1->data, 1, $v2->data, 1);
    }
    
    /**
     *  Copy the vector from X to Y.
     *    Y := X
     *  @param int $n           Number of elements in each vector
     *  @param \numphp\vector $vect_X        Vector X buffer
     *  @param \numphp\vector $vect_Y        Vector Y buffer
     *  @return void
     */
    public static function scopy(\numphp\vector $vect_X, \numphp\vector $vect_Y) {
        self::init();
        return self::$ffi_blas->cblas_scopy($vect_X->col, $vect_X->data, 1, 
                 $vect_Y->data, 1);
    }
    
    /**
     * Compute the Euclidean norm of a vector.
     *    ret := ||X||
     * @param \numphp\vector $v
     * @return float
     */
    public static function snrm2(\numphp\vector $v):float {
        self::init();
        return self::$ffi_blas->cblas_snrme($v->col, $v->data, 1);
    }
    /**
     *  Add vectors
     *    Y := alpha * X + Y
     *  @param int $n           Number of elements in each vector
     *  @param float $alpha     Coefficient of scalar multiple of X vector
     *  @param \numphp\vector $vect_X        Vector X buffer
     *  @param \numphp\vector $vect_Y        Vector Y buffer
     *  @return void
     */
    public static function axpy(float $alpha, \numphp\vector $vect_X, \numphp\vector $vect_Y) {
        self::init();
        return self::$ffi_blas->cblas_saxpy($vect_X->col, $alpha, $vect_X->data,
                 1, $vect_Y->data, 1);
    }
    
    /**
     *  Calculates the sum of the absolute values of each component of the vector.
     *    ret := |x_1| + ... + |x_n|
     *  @param int $n           Number of elements in each vector
     *  @param \numphp\vector $vector        Vector X buffer
     *  @return float
     */
    public static function asum(\numphp\vector $vector): float {
        self::init();
        return self::$ffi_blas->cblas_sasum($vector->col, $vector->data, 1);
    }
    
    /**
     * Rotate about a given point.
     *    X(i) := c * X(i) + s * Y(i)
     *    Y(i) :=-s * X(i) + c * Y(i)
     * @param \numphp\vector $v1 Vector X buffer
     * @param \numphp\vector $v2 Vector Y buffer
     * @param float $c         value of cos A(Value obtained with rotg function.)
     * @param float $s         value of sin A(Value obtained with rotg function.)
     * 
     */
    public static function srot(\numphp\vector $v1, \numphp\vector $v2,float $c, float $s) {
        self::init();
        return self::$ffi_blas->cblas_srot($v1->col, $v1->data, 1,
                $v2->data, 1, $c, $s);
    }
    
    /**
     *  Give the point P (a, b).
     *  Rotate this point to givens and calculate the parameters a, b, c,
     *  and s to make the y coordinate zero.
     *    Conditions description:
     *       c * a + s * b = r
     *       -s * a + c * b = 0
     *       r = ||(a,b)||
     *       c^2 + s^2 = 1
     *       z=s if |a| > |b|
     *       z=1/c if |a| <= |b| and c != 0 and r != 0
     *    Find r, z, c, s that satisfies the above description.
     *    However, when r = 0, z = 0, c = 1, and s = 0 are returned.
     *    Also, if c = 0, | a | <= | b | and c! = 0 and r! = 0, z = 1 is returned.
     *  @param float $a     X-coordinate of P: The calculated r value is stored and returned
     *  @param float $b     Y-coordinate of P: The calculated z value is stored and returned
     *  @param float $c     Stores the calculated value of c
     *  @param float $s     Stores the calculated value of s
     *  @return void
     */
    
    public static function srotg(float $a, float $b, float $c, float $s) {
        self::init();
        return self::$ffi_blas->cblas_srotg($a, $b, $c, $s);
    }
    
    /**
     * Multiply vector by scalar.
     *    X := alpha * X
     * @param float $alpha     Coefficient of scalar multiple of V vector
     * @param \numphp\vector $v
     * 
     */
    public static function sscal(float $alpha, \numphp\vector $v) {
        self::init();
        return self::$ffi_blas->cblas_sscal($v->col, $alpha, $v->data, 1);
    }
    

}
