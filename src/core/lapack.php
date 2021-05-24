<?php
namespace numphp\core;

/**
 * php interface for LAPACK
 * 
 * @package NumPhp\Lapack
 * @category Scientific Computing
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */

class lapack {
    const ROW_MAJOR = 101, COL_MAJOR = 102;
    const Upper = 'U', Lower = 'L';
    public static $ffi_lapack;
    
    public static function init() {
        if (is_null(self::$ffi_lapack)) {
            self::$ffi_lapack = \FFI::load(__DIR__ . '/lapack.h');
        }
        return self::$ffi_lapack;
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $ipiv
     * @param int $matLayout
     * @return int
     */
    public static function sgetrf(\numphp\matrix $mat, \numphp\vector $ipiv ,int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_sgetrf($matLayout, $mat->row, $mat->col, $mat->data, $mat->row, $ipiv->data);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $ipiv
     * @param int $matLayout
     * @return int
     */
    public static function dgetrf(\numphp\matrix $mat, \numphp\vector $ipiv ,int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dgetrf($matLayout, $mat->row, $mat->col, $mat->data, $mat->row, $ipiv->data);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $ipiv
     * @param int $matLayout
     * @return int
     */
    public static function sgetri(\numphp\matrix $mat, \numphp\vector $ipiv, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_sgetri($matLayout, $mat->row, $mat->data, $mat->row, $ipiv->data);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $ipiv
     * @param int $matLayout
     * @return int
     */
    public static function dgetri(\numphp\matrix $mat, \numphp\vector $ipiv, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dgetri($matLayout, $mat->row, $mat->data, $mat->row, $ipiv->data);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $s
     * @param \numphp\matrix $u
     * @param \numphp\matrix $vt
     * @param int $matLayout
     * @return int
     */
    public static function sgesdd(\numphp\matrix $mat, \numphp\vector $s, \numphp\matrix $u, \numphp\matrix $vt, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_sgesdd($matLayout,'A',$mat->row,$mat->col,$mat->data,$mat->col,$s->data,$u->data,$mat->row,$vt->data,$mat->col);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $s
     * @param \numphp\matrix $u
     * @param \numphp\matrix $vt
     * @param int $matLayout
     * @return int
     */
    public static function dgesdd(\numphp\matrix $mat, \numphp\vector $s, \numphp\matrix $u, \numphp\matrix $vt, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dgesdd($matLayout,'A',$mat->row,$mat->col,$mat->data,$mat->col,$s->data,$u->data,$mat->row,$vt->data,$mat->col);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param string $uplo
     * @param int $matLayout
     * @return int
     */
    public static function spotrf(\numphp\matrix $mat, $uplo= self::Lower, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_spotrf($matLayout, $uplo, $mat->col, $mat->data, $mat->col);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param string $uplo
     * @param int $matLayout
     * @return int
     */
    public static function dpotrf(\numphp\matrix $mat, $uplo= self::Lower, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dpotrf($matLayout, $uplo, $mat->col, $mat->data, $mat->col);
    }

    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $wr
     * @param \numphp\vector $wi
     * @param \numphp\matrix $vr
     * @param int $matLayout
     * @return int
     */
    public static function sgeev(\numphp\matrix $mat, \numphp\vector $wr, \numphp\vector $wi, \numphp\matrix $vr, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_sgeev($matLayout,'N','V',$mat->col,$mat->data,$mat->col,$wr->data,$wi->data,null,$mat->col,$vr->data, $mat->col);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $wr
     * @param \numphp\vector $wi
     * @param \numphp\matrix $vr
     * @param int $matLayout
     * @return int
     */
    public static function dgeev(\numphp\matrix $mat, \numphp\vector $wr, \numphp\vector $wi, \numphp\matrix $vr, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dgeev($matLayout,'N','V',$mat->col,$mat->data,$mat->col,$wr->data,$wi->data,null,$mat->col,$vr->data, $mat->col);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $wr
     * @param int $matLayout
     * @return int
     */
    public static function ssyev(\numphp\matrix $mat, \numphp\vector $wr, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_ssyev($matLayout,'V', 'U', $mat->col, $mat->data, $mat->col, $wr->data);
    }
    
    /**
     * 
     * @param \numphp\matrix $mat
     * @param \numphp\vector $wr
     * @param int $matLayout
     * @return int
     */
    public static function dsyev(\numphp\matrix $mat, \numphp\vector $wr, int $matLayout = self::ROW_MAJOR) {
        self::init();
        return self::$ffi_lapack->LAPACKE_dsyev($matLayout,'V', 'U', $mat->col, $mat->data, $mat->col, $wr->data);
    }
    
    /**
     * 
     * @param string $norm
     * @param \numphp\matrix $m
     * @param int $matLayout
     * @return type
     */
    public static function slange(string $norm, \numphp\matrix $m, int $matLayout = self::ROW_MAJOR){
        self::init();
        return self::$ffi_lapack->LAPACKE_slange($matLayout, $norm, $m->row, $m->col, $m->data,$m->col);
    }
    
    /**
     * 
     * @param string $norm
     * @param \numphp\matrix $m
     * @param int $matLayout
     * @return type
     */
    public static function dlange(string $norm, \numphp\matrix $m, int $matLayout = self::ROW_MAJOR){
        self::init();
        return self::$ffi_lapack->LAPACKE_dlange($matLayout, $norm, $m->row, $m->col, $m->data,$m->col);
    }
    
    


    public static function sgels() {
        self::init();
        
    }
}
