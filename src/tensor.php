<?php

declare(strict_types=1);

namespace numphp;

use FFI,
    InvalidArgumentException;

if (!extension_loaded('FFI')) {
    die('FFI extension required');
}

/** A fast lite memory effiecent Scientific Computing in php
 * 
 * @category Scientific Computing
 * @package numphp/tensor
 * @author ghostjat
 * @version 0.0.1.alpha
 */
class tensor {

    const TWO_PI = 2. * M_PI, EPSILON = 1e-8;
    const CblasLeft = 141, CblasRight = 142;
    const CblasUpper = 121, CblasLower = 122;
    const CblasNonUnit = 131, CblasUnit = 132;
    const CblasRowMajor = 101, CblasColMajor = 102;
    const CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113;
    const Float = 1, Double = 2;

    public $tp, $row, $col, $view, $type;
    public static $ffi_blas = null, $_time = null, $_mem = null;

    protected function __construct(int $row, int $col=null) {
        if ($row < 1) {
            throw new InvalidArgumentException(' * To create Numphp/tensor row must be greater than 0!, Op Failed! * ');
        }
        if(is_null($col)) {
            $this->type = 1;
            $this->row = $row;
            $this->tp = self::c_Vector($this->row);
        }
        else {
            $this->type = 0;
            $this->row = $row;
            $this->col = $col;
            $this->tp = self::c_Matrix($this->row, $this->col);
        }
        return $this;
    }
    
    /**
     * create empty vector/2d matrix
     * @param int $row
     * @param int $col
     * @return \numphp\tensor
     */
    public static function factory(int $row,int $col = null): tensor {
        return new self($row,$col);
    }
    
    /**
     * create vector/2d matrix using php array
     * @param array $data
     * @return \numphp\tensor
     */
    public static function ar(array $data): tensor {
        if (is_array($data) && is_array($data[0])) {
            $ar = self::factory(count($data), count($data[0]));
            $ar->setData($data);
        } else {
            $row = count($data);
            $ar = self::factory($row, null);
            for ($i =0; $i < $row; ++$i) {
                $ar->tp[$i] = $data[$i];
            }
        }
        return $ar;
    }
    
    /**
     * Return Tensor with random values
     * @param int $row
     * @param int $col
     * @param const $type
     * @return \blas\tensor
     */
    public static function randn(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = rand() / $max;
            }
        }
        return $ar;
    }
    /**
     * Return vector with random values
     * @param int $row
     * @return \numphp\tensor
     */
    public static function randnVector(int $row) : tensor {
        $ar = self::factory($row);
        $max = getrandmax();
        for($i = 0; $i < $ar->row; ++$i) {
            $ar->tp[$i] = rand()/$max;
        }
        return $ar;
    }
    
    /**
     * create a 2d uniform matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function uniform(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = rand(-$max, $max) / $max;
            }
        }
        return $ar;
    }
    
    /**
     * cretae zero like 2d matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function zero(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = 0;
            }
        }
        return $ar;
    }
    
    /**
     * 
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function ones(int $row, int $col): tensor {
        $ar = self::factory($row, $col, self::Double);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = 1;
            }
        }
        return $ar;
    }

    /**
     * 
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function null(int $row, int $col): tensor {
        $ar = self::factory($row, $col, self::Double);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = null;
            }
        }
        return $ar;
    }
    
    
}
